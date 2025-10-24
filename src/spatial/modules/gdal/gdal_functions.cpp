#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/function/copy_function.hpp"


#include "cpl_string.h"
#include "cpl_vsi.h"
#include "cpl_vsi_error.h"
#include "cpl_vsi_virtual.h"

#include "gdal.h"
#include "ogr_core.h"
#include "ogr_api.h"

#include "duckdb/common/arrow/arrow_converter.hpp"
#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {
namespace {

//======================================================================================================================
// GDAL READ
//======================================================================================================================
namespace gdal_read {
//----------------------------------------------------------------------------------------------------------------------
// BIND
//----------------------------------------------------------------------------------------------------------------------
class BindData final : public TableFunctionData {
public:
	string file_path;
};

auto Bind(ClientContext &ctx, TableFunctionBindInput &input, vector<LogicalType> &col_types, vector<string> &col_names)
	-> unique_ptr<FunctionData> {

	auto result = make_uniq<BindData>();

	result->file_path = input.inputs[0].GetValue<string>();

	const auto dataset = GDALOpenEx(result->file_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, nullptr, nullptr, nullptr);
	if (!dataset) {
		GDALClose(dataset);
		throw IOException("Could not open GDAL dataset at: %s", result->file_path);
	}

	if (GDALDatasetGetLayerCount(dataset) <= 0) {
		GDALClose(dataset);
		throw IOException("GDAL dataset contains no layers at: %s", result->file_path);
	}

	const auto layer = GDALDatasetGetLayer(dataset, 0);
	if (!layer) {
		GDALClose(dataset);
		throw IOException("Could not get GDAL layer at: %s", result->file_path);
	}

	ArrowArrayStream stream;
	if (!OGR_L_GetArrowStream(layer, &stream, nullptr)) {
		GDALClose(dataset);
		throw IOException("Could not get GDAL Arrow stream at: %s", result->file_path);
	}

	ArrowSchema schema;
	if (stream.get_schema(&stream, &schema) != 0) {
		stream.release(&stream);
		GDALClose(dataset);
		throw IOException("Could not get GDAL Arrow schema at: %s", result->file_path);
	}

	// Convert Arrow schema to DuckDB types
	for (int64_t i = 0; i < schema.n_children; i++) {
		auto &child_schema = *schema.children[i];
		const auto type = ArrowType::GetTypeFromSchema(ctx.db->config, child_schema);
		col_names.push_back(child_schema.name);
		col_types.push_back(type->GetDuckType());
	}

	// Release stream, schema and dataset
	schema.release(&schema);
	stream.release(&stream);
	GDALClose(dataset);


	return std::move(result);
}

class GlobalState final : public GlobalTableFunctionState {
public:

	~GlobalState() override {
		if (dataset) {
			GDALClose(dataset);
			dataset = nullptr;
		}

		if (stream.release) {
			stream.release(&stream);
		}
	}

	GDALDatasetH dataset;
	OGRLayerH layer;
	ArrowArrayStream stream;
	vector<unique_ptr<ArrowType>> col_types;
};

auto InitGlobal(ClientContext &context, TableFunctionInitInput &input) -> unique_ptr<GlobalTableFunctionState> {
	auto &data = input.bind_data->Cast<BindData>();

	const auto dataset = GDALOpenEx(data.file_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, nullptr, nullptr, nullptr);
	if (!dataset) {
		throw IOException("Could not open GDAL dataset at: foo");
	}

	auto result = make_uniq<GlobalState>();
	result->dataset = dataset;

	// Get the first layer
	result->layer = GDALDatasetGetLayer(dataset, 0);


	string str = "MAX_FEATURES_IN_BATCH=2048";
	vector<char> buf;
	buf.insert(buf.end(), str.begin(), str.end());
	buf.push_back('\0');
	vector<char*> layer_options;
	layer_options.push_back(buf.data());
	layer_options.push_back(nullptr);

	// Open the Arrow stream
	if (!OGR_L_GetArrowStream(result->layer, &result->stream, layer_options.data())) {
		GDALClose(dataset);
		throw IOException("Could not get GDAL Arrow stream at: foo");
	}

	ArrowSchema schema;
	if (result->stream.get_schema(&result->stream, &schema) != 0) {
		result->stream.release(&result->stream);
		GDALClose(dataset);
		throw IOException("Could not get GDAL Arrow schema at: foo");
	}

	// Store the column types
	for (int64_t i = 0; i < schema.n_children; i++) {
		auto &child_schema = *schema.children[i];
		result->col_types.push_back(ArrowType::GetTypeFromSchema(context.db->config, child_schema));
	}

	return std::move(result);
}

void Scan(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &bdata = input.bind_data->Cast<BindData>();
	auto &state = input.global_state->Cast<GlobalState>();

	ArrowArray arrow_array;
	if (state.stream.get_next(&state.stream, &arrow_array) != 0 || arrow_array.release == nullptr) {
		// Finished reading
		output.SetCardinality(0);
		return;
	}

	// Now convert the Arrow array to DuckDB
	for (idx_t i = 0; i < arrow_array.n_children; i++) {
		auto &arr = *arrow_array.children[i];
		auto &vec = output.data[i];

		auto &arrow_type = *state.col_types[i];
		auto array_state = ArrowArrayScanState(context);
		// We need to make sure that our chunk will hold the ownership
		array_state.owned_data = duckdb::make_shared_ptr<duckdb::ArrowArrayWrapper>();
		array_state.owned_data->arrow_array = arrow_array;

		// We set it to nullptr to effectively transfer the ownership
		arrow_array.release = nullptr;

		switch (arrow_type.GetPhysicalType()) {
		case ArrowArrayPhysicalType::DICTIONARY_ENCODED:
			ArrowToDuckDBConversion::ColumnArrowToDuckDBDictionary(vec, arr, 0, array_state,
																		   arrow_array.length, arrow_type);
			break;
		case ArrowArrayPhysicalType::RUN_END_ENCODED:
			ArrowToDuckDBConversion::ColumnArrowToDuckDBRunEndEncoded(vec, arr, 0, array_state,
				arrow_array.length, arrow_type);
			break;
		case ArrowArrayPhysicalType::DEFAULT:
			ArrowToDuckDBConversion::SetValidityMask(vec, arr, 0,
				arrow_array.length, arrow_array.offset, -1);
			ArrowToDuckDBConversion::ColumnArrowToDuckDB(vec, arr, 0, array_state,
																 arrow_array.length, arrow_type);
			break;
		default:
			throw NotImplementedException("ArrowArrayPhysicalType not recognized");
		}
	}

	output.SetCardinality(arrow_array.length);
}

void Register(ExtensionLoader &loader) {
	TableFunction read_func("gdal_read", {LogicalType::VARCHAR}, Scan, Bind, InitGlobal);
	loader.RegisterFunction(read_func);
}

} // namespace gdal_read
//======================================================================================================================
// GDAL COPY
//======================================================================================================================

namespace gdal_copy {

//----------------------------------------------------------------------------------------------------------------------
// Bind
//----------------------------------------------------------------------------------------------------------------------
class BindData final : public TableFunctionData {
public:
	string file_path;
	string driver_name;
	string layer_name;
	vector<string> driver_options;
	vector<string> layer_options;
	string target_srs;
	OGRwkbGeometryType geometry_type;

	// Arrow info
	ClientProperties props;
	ArrowSchema schema;
	unordered_map<idx_t, const shared_ptr<ArrowTypeExtensionData>> extension_type_cast;

	~BindData() override {
		if (schema.release) {
			schema.release(&schema);
		}
	}
};

bool MatchOption(const char* name, const pair<string, vector<Value>> &option, bool list = false) {
	if (StringUtil::CIEquals(name, option.first)) {
		if (option.second.empty()) {
			throw BinderException("GDAL COPY option '%s' requires a value", name);
		}
		if (!list) {
			if (option.second.size() != 1) {
				throw BinderException("GDAL COPY option '%s' only accepts a single value", name);
			}
			if (option.second.back().type().id() != LogicalTypeId::VARCHAR) {
				throw BinderException("GDAL COPY option '%s' must be a string", name);
			}
		} else {
			for (auto &val : option.second) {
				if (val.type().id() != LogicalTypeId::VARCHAR) {
					throw BinderException("GDAL COPY option '%s' must be a list of strings", name);
				}
			}
		}
		return true;
	}
	return false;
}

auto Bind(ClientContext &context, CopyFunctionBindInput &input, const vector<string> &names,
	const vector<LogicalType> &sql_types) -> unique_ptr<FunctionData> {
	auto result = make_uniq<BindData>();

	// Set file path
	result->file_path = input.info.file_path;

	// Parse options
	for (auto &option : input.info.options) {

		if (MatchOption("DRIVER", option)) {
			result->driver_name = option.second.back().GetValue<string>();
			continue;
		}

		if (MatchOption("LAYER_NAME", option)) {
			result->layer_name = option.second.back().GetValue<string>();
			continue;
		}

		if (MatchOption("SRS", option) || MatchOption("CRS", option)) {
			result->target_srs = option.second.back().GetValue<string>();
			continue;
		}

		if (MatchOption("GEOMETRY_TYPE", option)) {
			auto type = option.second.back().GetValue<string>();
			if (StringUtil::CIEquals(type, "POINT")) {
				result->geometry_type = wkbPoint;
			} else if (StringUtil::CIEquals(type, "LINESTRING")) {
				result->geometry_type = wkbLineString;
			} else if (StringUtil::CIEquals(type, "POLYGON")) {
				result->geometry_type = wkbPolygon;
			} else if (StringUtil::CIEquals(type, "MULTIPOINT")) {
				result->geometry_type = wkbMultiPoint;
			} else if (StringUtil::CIEquals(type, "MULTILINESTRING")) {
				result->geometry_type = wkbMultiLineString;
			} else if (StringUtil::CIEquals(type, "MULTIPOLYGON")) {
				result->geometry_type = wkbMultiPolygon;
			} else if (StringUtil::CIEquals(type, "GEOMETRYCOLLECTION")) {
				result->geometry_type = wkbGeometryCollection;
			} else {
				throw BinderException("Unsupported GEOMETRY_TYPE: '%s'", type);
			}
			continue;
		}

		if (MatchOption("LAYER_CREATION_OPTIONS", option, true)) {
			for (auto &val : option.second) {
				result->layer_options.push_back(val.GetValue<string>());
			}
			continue;
		}

		if (MatchOption("DATASET_CREATION_OPTIONS", option, true)) {
			for (auto &val : option.second) {
				result->driver_options.push_back(val.GetValue<string>());
			}
			continue;
		}

		throw BinderException("Unknown GDAL COPY option: '%s'", option.first);
	}

	// Check that options are valid
	if (result->driver_name.empty()) {
		throw BinderException("GDAL COPY option 'DRIVER' is required");
	}

	if (result->layer_name.empty()) {
		auto &fs = FileSystem::GetFileSystem(context);
		result->layer_name = fs.ExtractBaseName(result->file_path);
	}

	// Check the driver
	const auto driver = GDALGetDriverByName(result->driver_name.c_str());
	if (!driver) {
		throw BinderException("Could not find GDAL driver: " + result->driver_name);
	}

	// Try to get the file extension from the driver
	const auto file_ext = GDALGetMetadataItem(driver, GDAL_DMD_EXTENSIONS, nullptr);
	if (file_ext) {
		input.file_extension = file_ext;
	} else {
		const auto file_exts = GDALGetMetadataItem(driver, GDAL_DMD_EXTENSIONS, nullptr);
		const auto exts = StringUtil::Split(file_exts, ' ');
		if (!exts.empty()) {
			input.file_extension = exts[0];
		}
	}

	// Driver-specific checks
	if (result->driver_name == "OpenFileGDB" && result->geometry_type == wkbUnknown) {
		throw BinderException("OpenFileGDB requires 'GEOMETRY_TYPE' parameter to be set when writing!");
	}

	// Setup arrow schema
	result->props = context.GetClientProperties();
	result->extension_type_cast = duckdb::ArrowTypeExtensionData::GetExtensionTypes(context, sql_types);
	ArrowConverter::ToArrowSchema(&result->schema, sql_types, names, result->props);

	return std::move(result);
}

//----------------------------------------------------------------------------------------------------------------------
// Global State
//----------------------------------------------------------------------------------------------------------------------
class GlobalState final : public GlobalFunctionData {
public:

	~GlobalState() override {
		if (dataset) {
			GDALClose(dataset);
			dataset = nullptr;
		}

		if (array.release) {
			array.release(&array);
			array.release = nullptr;
		}
	}

	void Open(const BindData &data) {

		const auto driver = GDALGetDriverByName(data.driver_name.c_str());
		if (!driver) {
			throw InvalidInputException("Could not find GDAL driver: " + data.driver_name);
		}

		// Make CPL list for driver options
		vector<const char*> cpl_driver_options;
		for (auto &option : data.driver_options) {
			cpl_driver_options.push_back(option.c_str());
		}
		cpl_driver_options.push_back(nullptr);

		// Create Dataset
		dataset = GDALCreate(driver, data.file_path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
		if (!dataset) {
			throw IOException("Could not create GDAL dataset at: " + data.file_path);
		}

		// Make CPL list for layer options
		vector<const char*> cpl_layer_options;
		for (auto &option : data.layer_options) {
			cpl_layer_options.push_back(option.c_str());
		}
		cpl_layer_options.push_back(nullptr);

		// Create Layer
		layer = GDALDatasetCreateLayer(dataset, data.driver_name.c_str(), nullptr, wkbUnknown, nullptr);
		if (!layer) {
			throw IOException("Could not create GDAL layer in dataset at: " + data.file_path);
		}

		// Create fields for all children
		auto geometry_field_count = 0;
		for (auto i = 0; i < data.schema.n_children; i++) {
			const auto child_schema = data.schema.children[i];

			// Check if this is a geometry field
			if (child_schema->metadata != nullptr) {
				// TODO: Look for arrow metadata!
				geometry_field_count++;
				if (geometry_field_count > 1) {
					throw NotImplementedException("Multiple geometry fields not supported yet");
				}
			} else {
				// Register normal attribute
				if (!OGR_L_CreateFieldFromArrowSchema(layer, child_schema, nullptr)) {
					throw IOException("Could not create field in GDAL layer for column: " + string(child_schema->name));
				}
			}
		}
	}
public:
	mutex lock;
	GDALDatasetH dataset;
	OGRLayerH layer;
	ArrowArray array;
};

auto InitGlobal(ClientContext &context, FunctionData &bdata, const string &path) -> unique_ptr<GlobalFunctionData> {
	auto &bind_data = bdata.Cast<BindData>();
	auto result = make_uniq<GlobalState>();

	result->Open(bind_data);

	return std::move(result);
}


//----------------------------------------------------------------------------------------------------------------------
// Local State
//----------------------------------------------------------------------------------------------------------------------
class LocalState final : public LocalFunctionData {
public:
	// No-op, we don't need any local state for now
};

auto InitLocal(ExecutionContext &context, FunctionData &bind_data) -> unique_ptr<LocalFunctionData> {
	auto result = make_uniq<LocalState>();
	return std::move(result);
}

//----------------------------------------------------------------------------------------------------------------------
// Sink
//----------------------------------------------------------------------------------------------------------------------
void Sink(ExecutionContext &context, FunctionData &bdata_p, GlobalFunctionData &gstate_p, LocalFunctionData &lstate_p,
		DataChunk &input) {

	const auto &bdata = bdata_p.Cast<BindData>();
	auto &gstate = gstate_p.Cast<GlobalState>();

	// Lock
	lock_guard<mutex> guard(gstate.lock);

	auto &arrow_array = gstate.array;
	auto &arrow_schema = bdata.schema;

	ArrowConverter::ToArrowArray(input, &arrow_array, bdata.props, bdata.extension_type_cast);
	OGR_L_WriteArrowBatch(gstate.layer, &arrow_schema, &arrow_array, nullptr);

	if (arrow_array.release) {
		arrow_array.release(&arrow_array);
		arrow_array.release = nullptr;
	}
}

//----------------------------------------------------------------------------------------------------------------------
// Combine
//----------------------------------------------------------------------------------------------------------------------
void Combine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate,
						LocalFunctionData &lstate) {

}

//----------------------------------------------------------------------------------------------------------------------
// Finalize
//----------------------------------------------------------------------------------------------------------------------
void Finalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate) {

}

//----------------------------------------------------------------------------------------------------------------------
// Register
//----------------------------------------------------------------------------------------------------------------------
void Register(ExtensionLoader &loader) {
	CopyFunction info("GDAL");

	info.copy_to_bind = Bind;
	info.copy_to_initialize_local = InitLocal;
	info.copy_to_initialize_global = InitGlobal;
	info.copy_to_sink = Sink;
	info.copy_to_combine = Combine;
	info.copy_to_finalize = Finalize;
	info.extension = "gdal";

	loader.RegisterFunction(info);
}

} // namespace gdal_copy
} // namespace

void RegisterExtraFunction(ExtensionLoader &loader) {
	gdal_copy::Register(loader);
	gdal_read::Register(loader);
}
} // namespace duckdb