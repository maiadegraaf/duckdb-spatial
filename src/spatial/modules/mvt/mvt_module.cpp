// Mapbox Vector Tiles (MVT) implementation

#include "spatial/modules/mvt/mvt_module.hpp"

#include "duckdb/common/types/hash.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "spatial/geometry/geometry_serialization.hpp"
#include "spatial/geometry/sgl.hpp"
#include "spatial/spatial_types.hpp"
#include "spatial/util/function_builder.hpp"

#include "protozero/buffer_vector.hpp"
#include "protozero/basic_pbf_writer.hpp"
#include "spatial/util/binary_reader.hpp"

namespace duckdb {

namespace {

//======================================================================================================================
// LocalState
//======================================================================================================================

class LocalState final : public FunctionLocalState {
public:
	explicit LocalState(ClientContext &context) : arena(BufferAllocator::Get(context)), allocator(arena) {
	}

	static unique_ptr<FunctionLocalState> Init(ExpressionState &state, const BoundFunctionExpression &expr,
	                                           FunctionData *bind_data);
	static LocalState &ResetAndGet(ExpressionState &state);

	string_t Serialize(Vector &vector, const sgl::geometry &geom);

	GeometryAllocator &GetAllocator() {
		return allocator;
	}

private:
	ArenaAllocator arena;
	GeometryAllocator allocator;
};

unique_ptr<FunctionLocalState> LocalState::Init(ExpressionState &state, const BoundFunctionExpression &expr,
                                                FunctionData *bind_data) {
	return make_uniq_base<FunctionLocalState, LocalState>(state.GetContext());
}

LocalState &LocalState::ResetAndGet(ExpressionState &state) {
	auto &local_state = ExecuteFunctionState::GetFunctionState(state)->Cast<LocalState>();
	local_state.arena.Reset();
	return local_state;
}

string_t LocalState::Serialize(Vector &vector, const sgl::geometry &geom) {
	const auto size = Serde::GetRequiredSize(geom);
	auto blob = StringVector::EmptyString(vector, size);
	Serde::Serialize(geom, blob.GetDataWriteable(), size);
	blob.Finalize();
	return blob;
}

//======================================================================================================================
// ST_TileEnvelope
//======================================================================================================================

struct ST_TileEnvelope {
	static constexpr double RADIUS = 6378137.0;
	static constexpr double PI = 3.141592653589793;
	static constexpr double CIRCUMFERENCE = 2 * PI * RADIUS;

	static void ExecuteWebMercator(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		TernaryExecutor::Execute<int32_t, int32_t, int32_t, string_t>(
		    args.data[0], args.data[1], args.data[2], result, args.size(),
		    [&](int32_t tile_zoom, int32_t tile_x, int32_t tile_y) {
			    validate_tile_zoom_argument(tile_zoom);
			    uint32_t zoom_extent = 1u << tile_zoom;
			    validate_tile_index_arguments(zoom_extent, tile_x, tile_y);
			    sgl::geometry bbox;
			    get_tile_bbox(lstate.GetAllocator(), zoom_extent, tile_x, tile_y, bbox);
			    return lstate.Serialize(result, bbox);
		    });
	}

	static void validate_tile_zoom_argument(int32_t tile_zoom) {
		if ((tile_zoom < 0) || (tile_zoom > 30)) {
			throw InvalidInputException("ST_TileEnvelope: tile_zoom must be in the range [0,30]");
		}
	}

	static void validate_tile_index_arguments(uint32_t zoom_extent, int32_t tile_x, int32_t tile_y) {
		if ((tile_x < 0) || (static_cast<uint32_t>(tile_x) >= zoom_extent)) {
			throw InvalidInputException("ST_TileEnvelope: tile_x is out of range for specified tile_zoom");
		}
		if ((tile_y < 0) || (static_cast<uint32_t>(tile_y) >= zoom_extent)) {
			throw InvalidInputException("ST_TileEnvelope: tile_y is out of range for specified tile_zoom");
		}
	}

	static void get_tile_bbox(GeometryAllocator &allocator, uint32_t zoom_extent, int32_t tile_x, int32_t tile_y,
	                          sgl::geometry &bbox) {
		double single_tile_width = CIRCUMFERENCE / zoom_extent;
		double single_tile_height = CIRCUMFERENCE / zoom_extent;
		double tile_left = get_tile_left(tile_x, single_tile_width);
		double tile_right = tile_left + single_tile_width;
		double tile_top = get_tile_top(tile_y, single_tile_height);
		double tile_bottom = tile_top - single_tile_height;

		sgl::polygon::init_from_bbox(allocator, tile_left, tile_bottom, tile_right, tile_top, bbox);
	}

	static double get_tile_left(uint32_t tile_x, double single_tile_width) {
		return -0.5 * CIRCUMFERENCE + (tile_x * single_tile_width);
	}

	static double get_tile_top(uint32_t tile_y, double single_tile_height) {
		return 0.5 * CIRCUMFERENCE - (tile_y * single_tile_height);
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
        The `ST_TileEnvelope` scalar function generates tile envelope rectangular polygons from specified zoom level and tile indices.

        This is used in MVT generation to select the features corresponding to the tile extent. The envelope is in the Web Mercator
        coordinate reference system (EPSG:3857). The tile pyramid starts at zoom level 0, corresponding to a single tile for the
        world. Each zoom level doubles the number of tiles in each direction, such that zoom level 1 is 2 tiles wide by 2 tiles high,
        zoom level 2 is 4 tiles wide by 4 tiles high, and so on. Tile indices start at `[x=0, y=0]` at the top left, and increase
        down and right. For example, at zoom level 2, the top right tile is `[x=3, y=0]`, the bottom left tile is `[x=0, y=3]`, and
        the bottom right is `[x=3, y=3]`.

        ```sql
        SELECT ST_TileEnvelope(2, 3, 1);
        ```
    )";
	static constexpr auto EXAMPLE = R"(
        SELECT ST_TileEnvelope(2, 3, 1);
        ┌───────────────────────────────────────────────────────────────────────────────────────────────────────────┐
        │                                         st_tileenvelope(2, 3, 1)                                          │
        │                                                 geometry                                                  │
        ├───────────────────────────────────────────────────────────────────────────────────────────────────────────┤
        │ POLYGON ((1.00188E+07 0, 1.00188E+07 1.00188E+07, 2.00375E+07 1.00188E+07, 2.00375E+07 0, 1.00188E+07 0)) │
        └───────────────────────────────────────────────────────────────────────────────────────────────────────────┘
    )";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(ExtensionLoader &loader) {
		FunctionBuilder::RegisterScalar(loader, "ST_TileEnvelope", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("tile_zoom", LogicalType::INTEGER);
				variant.AddParameter("tile_x", LogicalType::INTEGER);
				variant.AddParameter("tile_y", LogicalType::INTEGER);
				variant.SetReturnType(GeoTypes::GEOMETRY());
				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteWebMercator);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "conversion");
		});
	}
};

//======================================================================================================================
// ST_AsMVT
//======================================================================================================================
enum class MVTValueType : uint32_t {
	INT = 1,
	FLOAT = 2,
	STRING = 3,
	BOOL = 4,
};

struct MVTValue {
	MVTValueType type;
	uint32_t size;
	union {
		int64_t int_value;
		double double_value;
		const char *string_value;
		bool bool_value;
	};
};

struct MVTValueEq {
	bool operator()(const MVTValue &a, const MVTValue &b) const {
		if (a.type != b.type) {
			return false;
		}
		switch (a.type) {
		case MVTValueType::INT:
			return a.int_value == b.int_value;
		case MVTValueType::FLOAT:
			return a.double_value == b.double_value;
		case MVTValueType::STRING:
			return (a.size == b.size) && (strncmp(a.string_value, b.string_value, a.size) == 0);
		case MVTValueType::BOOL:
			return a.bool_value == b.bool_value;
		}
		return false; // Should not reach here
	}
};

struct MVTValueHash {
	size_t operator()(const MVTValue &val) const {
		// Use duckdb::Hash
		size_t h1 = duckdb::Hash(static_cast<uint32_t>(val.type));
		size_t h2 = 0;
		switch (val.type) {
		case MVTValueType::INT:
			h2 = duckdb::Hash(val.int_value);
			break;
		case MVTValueType::FLOAT:
			h2 = duckdb::Hash(val.double_value);
			break;
		case MVTValueType::STRING:
			h2 = duckdb::Hash(val.string_value, val.size);
			break;
		case MVTValueType::BOOL:
			h2 = duckdb::Hash(val.bool_value);
			break;
		}
		return h1 ^ (h2 << 1); // Combine the two hashes
	}
};

using MVTValueDictionary = unordered_map<MVTValue, uint32_t, MVTValueHash, MVTValueEq>;

struct MVTFeature {
	MVTFeature *next;
	uint32_t id;
	uint32_t type;
	uint32_t geom_array_size;
	uint32_t tags_array_size;
	uint32_t *geom_array_data;
	uint32_t *tags_array_keys;
	MVTValue *tags_array_vals;
};

struct MVTLayer {
	MVTFeature *features_head = nullptr;
	MVTFeature *features_tail = nullptr;

	void Absorb(MVTLayer &other) {
		// Append other's features to this layer
		if (other.features_head) {
			if (features_tail) {
				features_tail->next = other.features_head;
				features_tail = other.features_tail;
			} else {
				features_head = other.features_head;
				features_tail = other.features_tail;
			}
			other.features_head = nullptr;
			other.features_tail = nullptr;
		}
	}

	void Combine(ArenaAllocator &allocator, const MVTLayer &other) {
		// Copy the features from the other into this, but reference the same values
		auto other_feature = other.features_head;
		while (other_feature) {
			const auto new_feature_mem = allocator.AllocateAligned(sizeof(MVTFeature));
			const auto new_feature = new (new_feature_mem) MVTFeature();

			// Copy the feature data
			*new_feature = *other_feature;

			new_feature->next = nullptr;
			if (features_tail) {
				features_tail->next = new_feature;
				features_tail = new_feature;
			} else {
				features_head = new_feature;
				features_tail = new_feature;
			}

			other_feature = other_feature->next;
		}
	}

	// Write the layer to the buffer
	void Finalize(const uint32_t extent, const vector<string> &tag_names, const string &layer_name,
	              vector<char> &buffer, MVTValueDictionary &tag_dict) {

		protozero::basic_pbf_writer<std::vector<char>> tile_writer {buffer};
		protozero::basic_pbf_writer<std::vector<char>> layer_writer {tile_writer, 3}; // layers = 3

		// Add version
		layer_writer.add_uint32(15, 2);

		// Layer name = 1
		layer_writer.add_string(1, layer_name);

		// Add layer name
		//layer_writer.add_string(1, bdata.layer_name);

		uint64_t fid = 0;

		auto feature = features_head;
		while (feature) {

			protozero::basic_pbf_writer<std::vector<char>> feature_writer {layer_writer, 2}; // features = 2

			// Id = 1
			feature_writer.add_uint64(1, fid++);

			// Tags = 2
			{
				protozero::detail::packed_field_varint<std::vector<char>, uint32_t> tags_writer(feature_writer, 2);
				for (uint32_t tag_idx = 0; tag_idx < feature->tags_array_size; tag_idx++) {
					const auto &key_idx = feature->tags_array_keys[tag_idx];
					const auto &val = feature->tags_array_vals[tag_idx];

					// Try to find the value in the dictionary
					// If it exists, we use the existing index
					// If it does not exist, we add it to the dictionary and use the newly added index
					const auto val_idx =
					    tag_dict.insert(make_pair(val, static_cast<uint32_t>(tag_dict.size()))).first->second;

					tags_writer.add_element(key_idx);
					tags_writer.add_element(val_idx);
				}
			}

			// Type = 3
			feature_writer.add_uint32(3, feature->type);

			// Geometry = 4
			feature_writer.add_packed_uint32(4, feature->geom_array_data,
			                                 feature->geom_array_data + feature->geom_array_size);

			feature = feature->next;
		}

		// Tag Keys = 3
		for (auto &key : tag_names) {
			layer_writer.add_string(3, key);
		}

		for (const auto &tag : tag_dict) {
			auto &val = tag.first;
			protozero::basic_pbf_writer<std::vector<char>> val_writer {layer_writer, 4}; // values = 4
			switch (val.type) {
			case MVTValueType::INT: {
				val_writer.add_int64(4, val.int_value);
			} break;
			case MVTValueType::FLOAT: {
				layer_writer.add_double(3, val.double_value);
			} break;
			case MVTValueType::STRING: {
				layer_writer.add_string(1, val.string_value, val.size);
			} break;
			default:
				throw InternalException("ST_AsMVT: Unsupported MVT value type");
			}
		}

		// Extent = 5
		layer_writer.add_uint32(5, extent);
	}
};

class MVTFeatureBuilder {
public:
	void Reset() {
		id = 0;
		geometry_type = 0;
		geometry.clear();
		tags.clear();
	}

	void SetGeometry(const string_t &geom_blob) {

		BinaryReader cursor(geom_blob.GetData(), geom_blob.GetSize());
		const auto type = static_cast<sgl::geometry_type>(cursor.Read<uint8_t>() + 1);
		const auto flags = cursor.Read<uint8_t>();
		cursor.Skip(sizeof(uint16_t));
		cursor.Skip(sizeof(uint32_t)); // padding

		// Parse flags
		const auto has_z = (flags & 0x01) != 0;
		const auto has_m = (flags & 0x02) != 0;
		const auto has_bbox = (flags & 0x04) != 0;

		const auto format_v1 = (flags & 0x40) != 0;
		const auto format_v0 = (flags & 0x80) != 0;

		if (format_v1 || format_v0) {
			// Unsupported version, throw an error
			throw NotImplementedException(
			    "This geometry seems to be written with a newer version of the DuckDB spatial library that is not "
			    "compatible with this version. Please upgrade your DuckDB installation.");
		}

		if (has_bbox) {
			// Skip past bbox if present
			cursor.Skip(sizeof(float) * 2 * (2 + has_z + has_m));
		}

		// Read the first type
		cursor.Skip(sizeof(uint32_t));

		const auto vertex_width = (2 + (has_z ? 1 : 0) + (has_m ? 1 : 0)) * sizeof(double);
		const auto vertex_space = vertex_width - (2 * sizeof(double)); // Space for x and y

		switch (type) {
		case sgl::geometry_type::POINT: {
			geometry_type = 1; // MVT_POINT

			// Read the point geometry
			const auto vertex_count = cursor.Read<uint32_t>();
			if (vertex_count == 0) {
				// No vertices, skip
				throw InvalidInputException("ST_AsMVT: POINT geometry cant be empty");
			}
			const auto x = CastDouble(cursor.Read<double>());
			const auto y = CastDouble(cursor.Read<double>());
			cursor.Skip(vertex_space); // Skip z and m if present

			geometry.push_back((1 & 0x7) | (1 << 3)); // MoveTo, 1 part
			geometry.push_back(protozero::encode_zigzag32(x));
			geometry.push_back(protozero::encode_zigzag32(y));

		} break;
		case sgl::geometry_type::LINESTRING: {
			geometry_type = 2; // MVT_LINESTRING

			const auto vertex_count = cursor.Read<uint32_t>();
			if (vertex_count < 2) {
				// Invalid linestring, skip
				throw InvalidInputException("ST_AsMVT: LINESTRING geometry cant contain less than 2 vertices");
			}
			// Read the vertices
			int32_t cursor_x = 0;
			int32_t cursor_y = 0;

			for (uint32_t vertex_idx = 0; vertex_idx < vertex_count; vertex_idx++) {

				const auto x = CastDouble(cursor.Read<double>());
				const auto y = CastDouble(cursor.Read<double>());
				cursor.Skip(vertex_space); // Skip z and m if present

				if (vertex_idx == 0) {
					geometry.push_back((1 & 0x7) | (1 << 3)); // MoveTo, 1 part
					geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
					geometry.push_back(protozero::encode_zigzag32(y - cursor_y));
					geometry.push_back((2 & 0x7) | ((vertex_count - 1) << 3)); // LineTo, part count
				} else {
					geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
					geometry.push_back(protozero::encode_zigzag32(y - cursor_y));
				}

				cursor_x = x;
				cursor_y = y;
			}
		} break;
		case sgl::geometry_type::POLYGON: {
			geometry_type = 3; // MVT_POLYGON

			const auto part_count = cursor.Read<uint32_t>();
			if (part_count == 0) {
				// No parts, invalid
				throw InvalidInputException("ST_AsMVT: POLYGON geometry cant be empty");
			}

			int32_t cursor_x = 0;
			int32_t cursor_y = 0;

			auto ring_cursor = cursor;
			cursor.Skip((part_count * 4) + (part_count % 2 == 1 ? 4 : 0)); // Skip part types and padding
			for (uint32_t part_idx = 0; part_idx < part_count; part_idx++) {
				const auto vertex_count = ring_cursor.Read<uint32_t>();
				if (vertex_count < 3) {
					// Invalid polygon, skip
					throw InvalidInputException("ST_AsMVT: POLYGON ring cant contain less than 3 vertices");
				}

				for (uint32_t vertex_idx = 0; vertex_idx < vertex_count; vertex_idx++) {
					const auto x = CastDouble(cursor.Read<double>());
					const auto y = CastDouble(cursor.Read<double>());
					cursor.Skip(vertex_space); // Skip z and m if present

					if (vertex_idx == 0) {
						geometry.push_back((1 & 0x7) | (1 << 3)); // MoveTo, 1 part
						geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
						geometry.push_back(protozero::encode_zigzag32(y - cursor_y));
						geometry.push_back((2 & 0x7) | ((vertex_count - 2) << 3));

						cursor_x = x;
						cursor_y = y;

					} else if (vertex_idx == vertex_count - 1) {
						// Close the ring
						geometry.push_back((7 & 0x7) | (1 << 3)); // ClosePath
					} else {
						// Add the vertex
						geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
						geometry.push_back(protozero::encode_zigzag32(y - cursor_y));

						cursor_x = x;
						cursor_y = y;
					}
				}
			}
		} break;
		case sgl::geometry_type::MULTI_POINT: {
			geometry_type = 1; // MVT_POINT

			const auto part_count = cursor.Read<uint32_t>();
			if (part_count == 0) {
				throw InvalidInputException("ST_AsMVT: MULTI_POINT geometry cant be empty");
			}

			int32_t cursor_x = 0;
			int32_t cursor_y = 0;

			geometry.push_back((1 & 0x7) | (part_count << 3)); // MoveTo, part count

			// Read the parts
			for (uint32_t part_idx = 0; part_idx < part_count; part_idx++) {
				cursor.Skip(sizeof(uint32_t)); // Skip part type
				const auto vertex_count = cursor.Read<uint32_t>();
				if (vertex_count == 0) {
					// No vertices, skip
					throw InvalidInputException("ST_AsMVT: POINT geometry cant be empty");
				}

				const auto x = CastDouble(cursor.Read<double>());
				const auto y = CastDouble(cursor.Read<double>());
				cursor.Skip(vertex_space); // Skip z and m if present

				geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
				geometry.push_back(protozero::encode_zigzag32(y - cursor_y));

				cursor_x = x;
				cursor_y = y;
			}
		} break;
		case sgl::geometry_type::MULTI_LINESTRING: {
			geometry_type = 2; // MVT_LINESTRING

			// Read the multi-linestring geometry
			const auto part_count = cursor.Read<uint32_t>();
			if (part_count == 0) {
				// No parts, invalid
				throw InvalidInputException("ST_AsMVT: MULTI_LINESTRING geometry cant be empty");
			}
			int32_t cursor_x = 0;
			int32_t cursor_y = 0;

			for (uint32_t part_idx = 0; part_idx < part_count; part_idx++) {
				cursor.Skip(sizeof(uint32_t)); // Skip part type
				const auto vertex_count = cursor.Read<uint32_t>();

				if (vertex_count < 2) {
					// Invalid linestring, skip
					throw InvalidInputException("ST_AsMVT: LINESTRING geometry cant contain less than 2 vertices");
				}

				for (uint32_t vertex_idx = 0; vertex_idx < vertex_count; vertex_idx++) {

					const auto x = CastDouble(cursor.Read<double>());
					const auto y = CastDouble(cursor.Read<double>());
					cursor.Skip(vertex_space); // Skip z and m if present

					if (vertex_idx == 0) {
						geometry.push_back((1 & 0x7) | (1 << 3)); // MoveTo, 1 part
						geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
						geometry.push_back(protozero::encode_zigzag32(y - cursor_y));
						geometry.push_back((2 & 0x7) | ((vertex_count - 2) << 3)); // LineTo, part count
					} else {
						geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
						geometry.push_back(protozero::encode_zigzag32(y - cursor_y));
					}

					cursor_x = x;
					cursor_y = y;
				}
			}
		} break;
		case sgl::geometry_type::MULTI_POLYGON: {
			geometry_type = 3; // MVT_POLYGON

			// Read the multi-linestring geometry
			const auto poly_count = cursor.Read<uint32_t>();
			if (poly_count == 0) {
				// No parts, invalid
				throw InvalidInputException("ST_AsMVT: MULTI_POLYGON geometry cant be empty");
			}

			int32_t cursor_x = 0;
			int32_t cursor_y = 0;

			for (uint32_t poly_idx = 0; poly_idx < poly_count; poly_idx++) {
				cursor.Skip(sizeof(uint32_t)); // Skip part type
				const auto part_count = cursor.Read<uint32_t>();
				if (part_count == 0) {
					// No parts, invalid
					throw InvalidInputException("ST_AsMVT: POLYGON geometry cant be empty");
				}

				auto ring_cursor = cursor;
				cursor.Skip((part_count * 4) + (part_count % 2 == 1 ? 4 : 0)); // Skip part types and padding

				for (uint32_t part_idx = 0; part_idx < part_count; part_idx++) {
					const auto vertex_count = ring_cursor.Read<uint32_t>();
					if (vertex_count < 3) {
						// Invalid polygon, skip
						throw InvalidInputException("ST_AsMVT: POLYGON ring cant contain less than 3 vertices");
					}

					for (uint32_t vertex_idx = 0; vertex_idx < vertex_count; vertex_idx++) {
						const auto x = CastDouble(cursor.Read<double>());
						const auto y = CastDouble(cursor.Read<double>());
						cursor.Skip(vertex_space); // Skip z and m if present

						if (vertex_idx == 0) {
							geometry.push_back((1 & 0x7) | (1 << 3)); // MoveTo, 1 part
							geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
							geometry.push_back(protozero::encode_zigzag32(y - cursor_y));
							geometry.push_back((2 & 0x7) | ((vertex_count - 2) << 3));

							cursor_x = x;
							cursor_y = y;

						} else if (vertex_idx == vertex_count - 1) {
							// Close the ring
							geometry.push_back((7 & 0x7) | (1 << 3)); // ClosePath
						} else {
							// Add the vertex
							geometry.push_back(protozero::encode_zigzag32(x - cursor_x));
							geometry.push_back(protozero::encode_zigzag32(y - cursor_y));

							cursor_x = x;
							cursor_y = y;
						}
					}
				}
			}
		} break;
		default:
			throw InvalidInputException("ST_AsMVT: unsupported geometry type %d", static_cast<int>(type));
		}
	}

	void AddProperty(uint32_t key, const string_t &value) {

		MVTValue v;
		v.type = MVTValueType::STRING;
		v.size = static_cast<uint32_t>(value.GetSize());
		v.string_value = value.GetData();

		tags.emplace_back(key, v);
	}

	void AddProperty(uint32_t key, int64_t value) {
		MVTValue v;
		v.type = MVTValueType::INT;
		v.size = sizeof(int64_t);
		v.int_value = value;

		tags.emplace_back(key, v);
	}

	bool IsEmpty() const {
		return geometry.empty();
	}

	void Finalize(ArenaAllocator &arena, MVTLayer &layer) {
		if (geometry.empty()) {
			// No geometry, skip
			return;
		}

		const auto feature_mem = arena.AllocateAligned(sizeof(MVTFeature));
		const auto feature_ptr = new (feature_mem) MVTFeature();

		feature_ptr->next = nullptr;
		feature_ptr->id = id;
		feature_ptr->type = geometry_type;

		// Copy over the geometry data
		feature_ptr->geom_array_data =
		    reinterpret_cast<uint32_t *>(arena.AllocateAligned(geometry.size() * sizeof(uint32_t)));
		feature_ptr->geom_array_size = static_cast<uint32_t>(geometry.size());
		memcpy(feature_ptr->geom_array_data, geometry.data(), geometry.size() * sizeof(uint32_t));

		// Copy over the tags
		feature_ptr->tags_array_size = static_cast<uint32_t>(tags.size());
		if (feature_ptr->tags_array_size != 0) {

			feature_ptr->tags_array_keys =
			    reinterpret_cast<uint32_t *>(arena.AllocateAligned(feature_ptr->tags_array_size * sizeof(uint32_t)));
			feature_ptr->tags_array_vals =
			    reinterpret_cast<MVTValue *>(arena.AllocateAligned(feature_ptr->tags_array_size * sizeof(MVTValue)));

			for (idx_t i = 0; i < tags.size(); i++) {
				feature_ptr->tags_array_keys[i] = tags[i].first;
				feature_ptr->tags_array_vals[i] = tags[i].second;
			}
		}

		// Append to the layer
		if (layer.features_tail) {
			layer.features_tail->next = feature_ptr;
			layer.features_tail = feature_ptr;
		} else {
			layer.features_head = feature_ptr;
			layer.features_tail = feature_ptr;
		}
	}

private:
	static int32_t CastDouble(double d) {
		if (d < static_cast<double>(std::numeric_limits<int32_t>::min()) ||
		    d > static_cast<double>(std::numeric_limits<int32_t>::max())) {
			throw InvalidInputException("ST_AsMVT: coordinate out of range for int32_t");
		}
		return static_cast<int32_t>(d);
	}

	uint32_t id = 0;
	uint32_t geometry_type = 0;
	vector<uint32_t> geometry;
	vector<pair<uint32_t, MVTValue>> tags;
};

struct ST_AsMVT {

	//------------------------------------------------------------------------------------------------------------------
	// Bind
	//------------------------------------------------------------------------------------------------------------------
	struct BindData final : FunctionData {

		idx_t geometry_column_idx = 0;
		string layer_name = "layer";
		uint32_t extent = 4096;
		vector<string> tag_names;

		unique_ptr<FunctionData> Copy() const override {
			auto result = make_uniq<BindData>();
			result->geometry_column_idx = geometry_column_idx;
			return std::move(result);
		}

		bool Equals(const FunctionData &other_p) const override {
			auto &other = other_p.Cast<BindData>();
			return geometry_column_idx == other.geometry_column_idx;
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, AggregateFunction &function,
	                                     vector<unique_ptr<Expression>> &arguments) {
		auto result = make_uniq<BindData>();

		string geom_name;

		// Figure part of the row is the geometry column
		const auto &row_type = arguments[0]->return_type;
		if (row_type.id() != LogicalTypeId::STRUCT) {
			throw InvalidInputException("ST_AsMVT: first argument must be a STRUCT (i.e. a row type)");
		}

		optional_idx geom_idx = optional_idx::Invalid();

		if (geom_name.empty()) {
			// Look for the first geometry column
			for (idx_t i = 0; i < StructType::GetChildCount(row_type); i++) {
				auto &child = StructType::GetChildType(row_type, i);
				if (child == GeoTypes::GEOMETRY()) {
					if (geom_idx != optional_idx::Invalid()) {
						throw InvalidInputException("ST_AsMVT: only one geometry column is allowed in the input row");
					}
					geom_idx = i;
				}
			}
		} else {
			// Look for the geometry column by name
			for (idx_t i = 0; i < StructType::GetChildCount(row_type); i++) {
				auto &child = StructType::GetChildType(row_type, i);
				auto &child_name = StructType::GetChildName(row_type, i);
				if (child == GeoTypes::GEOMETRY() && child_name == geom_name) {
					if (geom_idx != optional_idx::Invalid()) {
						throw InvalidInputException("ST_AsMVT: only one geometry column is allowed in the input row");
					}
					geom_idx = i;
				}
			}
		}
		if (!geom_idx.IsValid()) {
			throw InvalidInputException("ST_AsMVT: input row must contain a geometry column");
		}

		result->geometry_column_idx = geom_idx.GetIndex();

		return std::move(result);
	}

	//------------------------------------------------------------------------------------------------------------------
	// Initialize
	//------------------------------------------------------------------------------------------------------------------
	struct State {
		MVTLayer layer;
	};

	static idx_t StateSize(const AggregateFunction &) {
		return sizeof(State);
	}

	static void Initialize(const AggregateFunction &, data_ptr_t state_mem) {
		new (state_mem) State();
	}

	//------------------------------------------------------------------------------------------------------------------
	// Update
	//------------------------------------------------------------------------------------------------------------------
	static void Update(Vector inputs[], AggregateInputData &aggr, idx_t, Vector &state_vec, idx_t count) {
		const auto &bdata = aggr.bind_data->Cast<BindData>();
		const auto &row_cols = StructVector::GetEntries(inputs[0]);

		UnifiedVectorFormat state_format;
		UnifiedVectorFormat geom_format;
		vector<UnifiedVectorFormat> property_formats;
		vector<LogicalType> property_types;

		state_vec.ToUnifiedFormat(count, state_format);

		for (idx_t col_idx = 0; col_idx < row_cols.size(); col_idx++) {
			if (col_idx == bdata.geometry_column_idx) {
				row_cols[col_idx]->ToUnifiedFormat(count, geom_format);
			} else {
				property_formats.emplace_back();
				row_cols[col_idx]->ToUnifiedFormat(count, property_formats.back());
				property_types.push_back(row_cols[col_idx]->GetType());
			}
		}

		// Reusable geometry buffer
		MVTFeatureBuilder feature;

		for (idx_t row_idx = 0; row_idx < count; row_idx++) {
			const auto state_idx = state_format.sel->get_index(row_idx);
			auto &layer = UnifiedVectorFormat::GetData<State *>(state_format)[state_idx]->layer;

			const auto geom_idx = geom_format.sel->get_index(row_idx);
			if (!geom_format.validity.RowIsValid(geom_idx)) {
				// Skip if geometry is NULL
				continue;
			}

			auto &geom_blob = UnifiedVectorFormat::GetData<string_t>(geom_format)[geom_idx];

			// Reset the feature
			feature.Reset();

			// Set geometry
			feature.SetGeometry(geom_blob);

			// Add properties
			for (idx_t prop_vec_idx = 0; prop_vec_idx < property_formats.size(); prop_vec_idx++) {
				const auto &prop_format = property_formats[prop_vec_idx];
				const auto prop_row_idx = prop_format.sel->get_index(row_idx);
				if (!prop_format.validity.RowIsValid(prop_row_idx)) {
					// Skip if property is NULL
					continue;
				}

				// Switch on property type
				auto &prop_type = property_types[prop_vec_idx];
				switch (prop_type.id()) {
				case LogicalTypeId::VARCHAR: {
					auto &prop_val = UnifiedVectorFormat::GetData<string_t>(prop_format)[prop_row_idx];
					feature.AddProperty(static_cast<uint32_t>(prop_vec_idx), prop_val);
				} break;
				case LogicalTypeId::BIGINT: {
					auto &prop_val = UnifiedVectorFormat::GetData<int64_t>(prop_format)[prop_row_idx];
					feature.AddProperty(static_cast<uint32_t>(prop_vec_idx), prop_val);

				} break;
				default:
					throw InvalidInputException("ST_AsMVT: unsupported property type: %s", prop_type.ToString());
				}
			}

			feature.Finalize(aggr.allocator, layer);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Combine
	//------------------------------------------------------------------------------------------------------------------
	static void Combine(Vector &source_vec, Vector &target_vec, AggregateInputData &aggr, idx_t count) {
		UnifiedVectorFormat source_format;
		source_vec.ToUnifiedFormat(count, source_format);

		const auto source_ptr = UnifiedVectorFormat::GetData<State *>(source_format);
		const auto target_ptr = FlatVector::GetData<State *>(target_vec);

		for (idx_t row_idx = 0; row_idx < count; row_idx++) {
			auto &source = *source_ptr[source_format.sel->get_index(row_idx)];
			auto &target = *target_ptr[row_idx];

			if (aggr.combine_type == AggregateCombineType::ALLOW_DESTRUCTIVE) {
				// Absorb the feature data from source into target
				target.layer.Absorb(source.layer);
			} else {
				// Append the feature data from source to target
				target.layer.Combine(aggr.allocator, source.layer);
			}
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Finalize
	//------------------------------------------------------------------------------------------------------------------
	static void Finalize(Vector &state_vec, AggregateInputData &aggr, Vector &result, idx_t count, idx_t offset) {
		const auto &bdata = aggr.bind_data->Cast<BindData>();

		UnifiedVectorFormat state_format;
		state_vec.ToUnifiedFormat(count, state_format);
		const auto state_ptr = UnifiedVectorFormat::GetData<State *>(state_format);

		vector<char> buffer;
		MVTValueDictionary tag_dict;

		for (idx_t raw_idx = 0; raw_idx < count; raw_idx++) {
			auto &state = *state_ptr[state_format.sel->get_index(raw_idx)];
			const auto out_idx = raw_idx + offset;

			buffer.clear();
			tag_dict.clear();

			state.layer.Finalize(bdata.extent, bdata.tag_names, bdata.layer_name, buffer, tag_dict);

			// Now we have the layer buffer, we can write it to the result vector
			const auto result_data = FlatVector::GetData<string_t>(result);
			result_data[out_idx] = StringVector::AddStringOrBlob(result, buffer.data(), buffer.size());
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(ExtensionLoader &loader) {
		AggregateFunction agg({LogicalTypeId::ANY}, LogicalType::BLOB, StateSize, Initialize, Update, Combine, Finalize,
		                      nullptr, Bind);

		FunctionBuilder::RegisterAggregate(loader, "ST_AsMVT", [&](AggregateFunctionBuilder &func) {
			func.SetFunction(agg);
			func.SetDescription("Makes a vector tile from a set of geometries");

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

} // namespace
//======================================================================================================================
//  Register
//======================================================================================================================
void RegisterMapboxVectorTileModule(ExtensionLoader &loader) {
	ST_TileEnvelope::Register(loader);
	ST_AsMVT::Register(loader);
}

} // namespace duckdb
