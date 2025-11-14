#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/common/arrow/arrow_converter.hpp"
#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/parser/expression/constant_expression.hpp"
#include "duckdb/parser/expression/function_expression.hpp"

#include "gdal.h"
#include "ogr_core.h"
#include "ogr_api.h"
#include "ogr_srs_api.h"

#include "cpl_string.h"
#include "cpl_vsi.h"
#include "cpl_vsi_error.h"
#include "cpl_vsi_virtual.h"
#include "duckdb/parser/tableref/table_function_ref.hpp"

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

	int layer_idx = 0;
	bool keep_wkb = false;

	CPLStringList layer_options;
	CPLStringList dataset_options;
	CPLStringList dataset_sibling;
	CPLStringList dataset_drivers;

	int64_t estimated_cardinality = 0;
	unordered_set<idx_t> geometry_columns = {};

	bool can_filter = false;
	bool has_extent = false;
	bool has_filter = false;
	OGREnvelope layer_extent;
	OGREnvelope layer_filter;

	OGRwkbGeometryType layer_type = wkbUnknown;
};

auto Bind(ClientContext &ctx, TableFunctionBindInput &input, vector<LogicalType> &col_types, vector<string> &col_names)
    -> unique_ptr<FunctionData> {

	auto result = make_uniq<BindData>();

	// Pass file path
	result->file_path = input.inputs[0].GetValue<string>();

	// Parse options
	const auto dataset_options_param = input.named_parameters.find("open_options");
	if (dataset_options_param != input.named_parameters.end()) {
		for (auto &param : ListValue::GetChildren(dataset_options_param->second)) {
			result->dataset_options.AddString(StringValue::Get(param).c_str());
		}
	}

	const auto drivers_param = input.named_parameters.find("allowed_drivers");
	if (drivers_param != input.named_parameters.end()) {
		for (auto &param : ListValue::GetChildren(drivers_param->second)) {
			result->dataset_drivers.AddString(StringValue::Get(param).c_str());
		}
	}

	const auto siblings_params = input.named_parameters.find("sibling_files");
	if (siblings_params != input.named_parameters.end()) {
		for (auto &param : ListValue::GetChildren(siblings_params->second)) {
			result->dataset_sibling.AddString(StringValue::Get(param).c_str());
		}
	}

	const auto keep_wkb_param = input.named_parameters.find("keep_wkb");
	if (keep_wkb_param != input.named_parameters.end()) {
		result->keep_wkb = BooleanValue::Get(keep_wkb_param->second);
	}

	// Set additional default GDAL default options

	// This for OSM, but we don't know if we are reading OSM until we open the dataset, so just always set it for now.
	//result->dataset_options.AddString("INTERLEAVED_READING=YES");

	// This is so taht we dont have to deal with chunking ourselves, let GDAL do it for us
	result->layer_options.AddString(StringUtil::Format("MAX_FEATURES_IN_BATCH=%d", STANDARD_VECTOR_SIZE).c_str());

	// We always want GeoArrow geometry which DuckDB knows how to convert to GEOMETRY type, unless `keep_wkb` is set
	if (!result->keep_wkb) {
		result->layer_options.AddString("GEOMETRY_METADATA_ENCODING=GEOARROW");
	}

	// Open the dataset and get the Arrow schema
	const auto dataset =
	    GDALOpenEx(result->file_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
	    	result->dataset_drivers,
	    	result->dataset_options,
	    	result->dataset_sibling);

	if (!dataset) {
		throw IOException("Could not open GDAL dataset at: %s", result->file_path);
	}

	ArrowSchema schema;
	ArrowArrayStream stream;

	try {

		const auto layer_count = GDALDatasetGetLayerCount(dataset);
		if (layer_count <= 0) {
			throw IOException("GDAL dataset contains no layers at: %s", result->file_path);
		}

		// Find layer
		const auto layer_param = input.named_parameters.find("layer");

		if (layer_param != input.named_parameters.end()) {
			if (layer_param->second.type() == LogicalType::INTEGER) {
				// Find layer by index
				const auto layer_idx = IntegerValue::Get(layer_param->second);
				if (layer_idx < 0) {
					throw BinderException("Layer index must be positive");
				}
				if (layer_idx > layer_count) {
					throw BinderException(
					    StringUtil::Format("Layer index out of range (%s > %s)", layer_idx, layer_count));
				}
				result->layer_idx = layer_idx;
			} else if (layer_param->second.type() == LogicalType::VARCHAR) {
				// Find layer by name
				const auto &layer_name = StringValue::Get(layer_param->second);
				auto found = false;
				for (int i = 0; i < layer_count; i++) {
					const auto layer = GDALDatasetGetLayer(dataset, i);
					if (!layer) {
						continue;
					}
					if (OGR_L_GetName(layer) == layer_name) {
						result->layer_idx = i;
						found = true;
						break;
					}
				}
				if (!found) {
					throw BinderException("Could not find layer with name: %s", layer_name);
				}
			}
		}

		// Get the layer by index
		const auto layer = GDALDatasetGetLayer(dataset, result->layer_idx);
		if (!layer) {
			throw IOException("Could not get GDAL layer at: %s", result->file_path);
		}

		// Estimate cardinality
		result->estimated_cardinality = OGR_L_GetFeatureCount(layer, 0);

		// Get extent (Only if spatial filter is not pushed down!)
		if (OGR_L_GetExtent(layer, &result->layer_extent, 0) == OGRERR_NONE) {
			result->has_extent = true;
		}

		// Check if fast spatial filtering is available
		if (OGR_L_TestCapability(layer, OLCFastSpatialFilter)) {
			result->can_filter = true;
		}

		// Get the layer geometry type if available
		result->layer_type = OGR_L_GetGeomType(layer);

		// Get the arrow stream
		if (!OGR_L_GetArrowStream(layer, &stream, result->layer_options.List())) {
			throw IOException("Could not get GDAL Arrow stream at: %s", result->file_path);
		}

		// And the schema
		if (stream.get_schema(&stream, &schema) != 0) {
			throw IOException("Could not get GDAL Arrow schema at: %s", result->file_path);
		}

		// Convert Arrow schema to DuckDB types
		for (int64_t i = 0; i < schema.n_children; i++) {
			auto &child_schema = *schema.children[i];
			const auto gdal_type = ArrowType::GetTypeFromSchema(ctx.db->config, child_schema);
			auto duck_type = gdal_type->GetDuckType();

			// Track geometry columns to compute stats later
			if (duck_type.id() == LogicalTypeId::GEOMETRY) {
				result->geometry_columns.insert(i);
			}

			col_names.push_back(child_schema.name);
			col_types.push_back(std::move(duck_type));
		}

	} catch (...) {
		// Release stream, schema and dataset
		if (schema.release) {
			schema.release(&schema);
		}
		if (stream.release) {
			stream.release(&stream);
		}
		if (dataset) {
			GDALClose(dataset);
		}
		// Re-throw exception
		throw;
	}

	if (schema.release) {
		schema.release(&schema);
	}
	if (stream.release) {
		stream.release(&stream);
	}
	if (dataset) {
		GDALClose(dataset);
	}

	return std::move(result);
}

//----------------------------------------------------------------------------------------------------------------------
// FILTER (EXPRESSION) PUSHDOWN
//----------------------------------------------------------------------------------------------------------------------
auto Pushdown(ClientContext &context, LogicalGet &get, FunctionData *bind_data, vector<unique_ptr<Expression>> &filters)
    -> void {

	auto &bdata = bind_data->Cast<BindData>();

	if (!bdata.can_filter) {
		return;
	}

	if (bdata.geometry_columns.size() != 1) {
		return; // Only optimize if there is a single geometry column
	}

	optional_idx geom_filter_idx = optional_idx::Invalid();

	for (idx_t expr_idx = 0; expr_idx < filters.size(); expr_idx++) {
		const auto &expr = filters[expr_idx];

		if (expr->GetExpressionType() != ExpressionType::BOUND_FUNCTION) {
			continue;
		}
		if (expr->return_type != LogicalType::BOOLEAN) {
			continue;
		}
		const auto &func = expr->Cast<BoundFunctionExpression>();
		if (func.children.size() != 2) {
			continue;
		}

		if (func.children[0]->return_type.id() != LogicalTypeId::GEOMETRY ||
		    func.children[1]->return_type.id() != LogicalTypeId::GEOMETRY) {
			continue;
		}

		// The set of geometry predicates that can be optimized using the bounding box
		static constexpr const char *geometry_predicates[2] = {"&&", "st_intersects_extent"};

		auto found = false;
		for (const auto &name : geometry_predicates) {
			if (StringUtil::CIEquals(func.function.name.c_str(), name)) {
				found = true;
				break;
			}
		}
		if (!found) {
			// Not a geometry predicate we can optimize
			continue;
		}

		const auto lhs_kind = func.children[0]->GetExpressionType();
		const auto rhs_kind = func.children[1]->GetExpressionType();

		const auto lhs_is_const =
		    lhs_kind == ExpressionType::VALUE_CONSTANT && rhs_kind == ExpressionType::BOUND_COLUMN_REF;
		const auto rhs_is_const =
		    rhs_kind == ExpressionType::VALUE_CONSTANT && lhs_kind == ExpressionType::BOUND_COLUMN_REF;

		if (lhs_is_const == rhs_is_const) {
			// Both sides are constant or both sides are column refs
			continue;
		}

		auto &constant_expr = func.children[lhs_is_const ? 0 : 1]->Cast<BoundConstantExpression>();
		auto &geometry_expr = func.children[lhs_is_const ? 1 : 0]->Cast<BoundColumnRefExpression>();

		if (constant_expr.value.type().id() != LogicalTypeId::GEOMETRY) {
			// Constant is not geometry
			continue;
		}
		if (constant_expr.value.IsNull()) {
			// Constant is NULL
			continue;
		}
		if (geometry_expr.alias != "geom") {
			// Not the geometry column
			continue;
		}

		auto geom_extent = GeometryExtent::Empty();
		auto geom_binary = string_t(StringValue::Get(constant_expr.value));

		if (Geometry::GetExtent(geom_binary, geom_extent)) {
			bdata.has_filter = true;
			bdata.layer_filter.MinX = geom_extent.x_min;
			bdata.layer_filter.MinY = geom_extent.y_min;
			bdata.layer_filter.MaxX = geom_extent.x_max;
			bdata.layer_filter.MaxY = geom_extent.y_max;
		}

		// Set the index so we can remove it later
		// We can __ONLY__ do this if the filter predicate is "&&" or "st_intersects_extent"
		// as other predicates may require exact geometry evaluation, the filter cannot be fully removed
		geom_filter_idx = expr_idx;
		break;
	}

	if (geom_filter_idx != optional_idx::Invalid()) {
		// Remove the filter from the list
		filters.erase_at(geom_filter_idx.GetIndex());
	}
}

//----------------------------------------------------------------------------------------------------------------------
// GLOBAL STATE
//----------------------------------------------------------------------------------------------------------------------
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
	CPLStringList layer_options;
	OGRLayerH layer;
	ArrowArrayStream stream;
	vector<unique_ptr<ArrowType>> col_types;
	atomic<idx_t> features_read = {0};
};

auto InitGlobal(ClientContext &context, TableFunctionInitInput &input) -> unique_ptr<GlobalTableFunctionState> {
	auto &bdata = input.bind_data->Cast<BindData>();

	const auto dataset =
	    GDALOpenEx(bdata.file_path.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
	    	bdata.dataset_drivers,
	    	bdata.dataset_options,
	    	bdata.dataset_sibling);

	if (!dataset) {
		throw IOException("Could not open GDAL dataset at: foo");
	}

	auto result = make_uniq<GlobalState>();
	result->dataset = dataset;
	result->layer_options = bdata.layer_options;

	const auto driver = GDALGetDatasetDriver(dataset);
	if (strcmp(GDALGetDriverShortName(driver), "OSM") != 0) {
		// Get the layer by index
		result->layer = GDALDatasetGetLayer(dataset, bdata.layer_idx);
	} else {
		// Special case for OSM, which requires sequential reading of layers
		const auto layer_count = GDALDatasetGetLayerCount(dataset);
		for (int i = 0; i < layer_count; i++) {
			result->layer = GDALDatasetGetLayer(dataset, i);
			if (i == bdata.layer_idx) {
				// desired layer found
				break;
			}

			// else scan through and empty the layer
			OGRFeatureH feature;
			while ((feature = OGR_L_GetNextFeature(result->layer)) != nullptr) {
				OGR_F_Destroy(feature);
			}
		}
	}

	// Set the filter, if we got one
	if (bdata.has_filter) {
		OGR_L_SetSpatialFilterRect(result->layer, bdata.layer_filter.MinX, bdata.layer_filter.MinY,
		                           bdata.layer_filter.MaxX, bdata.layer_filter.MaxY);
	}

	CPLStringList layer_options;
	layer_options.AddString(StringUtil::Format("MAX_FEATURES_IN_BATCH=%d", STANDARD_VECTOR_SIZE).data());
	layer_options.AddString("GEOMETRY_METADATA_ENCODING=GEOARROW");

	// Open the Arrow stream
	if (!OGR_L_GetArrowStream(result->layer, &result->stream, result->layer_options.List())) {
		GDALClose(dataset);
		throw IOException("Could not get GDAL Arrow stream");
	}

	ArrowSchema schema;
	if (result->stream.get_schema(&result->stream, &schema) != 0) {
		result->stream.release(&result->stream);
		GDALClose(dataset);
		throw IOException("Could not get GDAL Arrow schema");
	}

	// Store the column types
	for (int64_t i = 0; i < schema.n_children; i++) {
		auto &child_schema = *schema.children[i];
		result->col_types.push_back(ArrowType::GetTypeFromSchema(context.db->config, child_schema));
	}

	return std::move(result);
}

//----------------------------------------------------------------------------------------------------------------------
// SCAN
//----------------------------------------------------------------------------------------------------------------------
void Scan(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
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
		array_state.owned_data = make_shared_ptr<ArrowArrayWrapper>();
		array_state.owned_data->arrow_array = arrow_array;

		// We set it to nullptr to effectively transfer the ownership
		arrow_array.release = nullptr;

		switch (arrow_type.GetPhysicalType()) {
		case ArrowArrayPhysicalType::DICTIONARY_ENCODED:
			ArrowToDuckDBConversion::ColumnArrowToDuckDBDictionary(vec, arr, 0, array_state, arrow_array.length,
			                                                       arrow_type);
			break;
		case ArrowArrayPhysicalType::RUN_END_ENCODED:
			ArrowToDuckDBConversion::ColumnArrowToDuckDBRunEndEncoded(vec, arr, 0, array_state, arrow_array.length,
			                                                          arrow_type);
			break;
		case ArrowArrayPhysicalType::DEFAULT:
			ArrowToDuckDBConversion::SetValidityMask(vec, arr, 0, arrow_array.length, arrow_array.offset, -1);
			ArrowToDuckDBConversion::ColumnArrowToDuckDB(vec, arr, 0, array_state, arrow_array.length, arrow_type);
			break;
		default:
			throw NotImplementedException("ArrowArrayPhysicalType not recognized");
		}
	}

	state.features_read += arrow_array.length;
	output.SetCardinality(arrow_array.length);
}

//------------------------------------------------------------------------------------------------------------------
// CARDINALITY
//------------------------------------------------------------------------------------------------------------------
auto Cardinality(ClientContext &context, const FunctionData *data) -> unique_ptr<NodeStatistics> {
	auto &bdata = data->Cast<BindData>();
	auto result = make_uniq<NodeStatistics>();

	if (bdata.estimated_cardinality > -1) {
		result->has_estimated_cardinality = true;
		result->estimated_cardinality = bdata.estimated_cardinality;
		result->has_max_cardinality = true;
		result->max_cardinality = bdata.estimated_cardinality;
	}

	return result;
}

//----------------------------------------------------------------------------------------------------------------------
// STATISTICS
//----------------------------------------------------------------------------------------------------------------------
auto Statistics(ClientContext &context, const FunctionData *bind_data, column_t column_index)
    -> unique_ptr<BaseStatistics> {

	auto &bdata = bind_data->Cast<BindData>();

	// If we have an extent, and the column is a geometry column, we can provide min/max stats
	if (bdata.has_extent) {

		// Check if this is the only geometry column
		const auto is_geom_col = bdata.geometry_columns.find(column_index) != bdata.geometry_columns.end();
		const auto is_only_one = bdata.geometry_columns.size() == 1;
		const auto has_stats = bdata.has_extent || bdata.layer_type != wkbUnknown;

		if (is_geom_col && is_only_one && has_stats) {
			auto stats = GeometryStats::CreateUnknown(LogicalType::GEOMETRY());

			if (bdata.has_extent) {
				auto &extent = GeometryStats::GetExtent(stats);
				extent.x_min = bdata.layer_extent.MinX;
				extent.x_max = bdata.layer_extent.MaxX;
				extent.y_min = bdata.layer_extent.MinY;
				extent.y_max = bdata.layer_extent.MaxY;
			}

			const auto geom_type = bdata.layer_type % 1000;
			const auto vert_type = bdata.layer_type / 1000;

			if ((geom_type >= 1) && (geom_type <= 7) && (vert_type >= 0) && (vert_type <= 3)) {
				auto &types = GeometryStats::GetTypes(stats);
				types.Clear();
				types.AddWKBType(static_cast<int32_t>(geom_type));
			}

			return stats.ToUnique();
		}
	}

	return nullptr;
}

//----------------------------------------------------------------------------------------------------------------------
// PROGRESS
//----------------------------------------------------------------------------------------------------------------------
auto Progress(ClientContext &context, const FunctionData *b_data, const GlobalTableFunctionState *g_state) -> double {
	auto &bdata = b_data->Cast<BindData>();
	auto &gstate = g_state->Cast<GlobalState>();

	if (bdata.estimated_cardinality < 0) {
		return 0.0;
	}

	const auto count = static_cast<double>(gstate.features_read.load());
	const auto total = static_cast<double>(bdata.estimated_cardinality);

	return MinValue(100.0 * (total / count), 100.0);
}

//------------------------------------------------------------------------------------------------------------------
// REPLACEMENT SCAN
//------------------------------------------------------------------------------------------------------------------
auto ReplacementScan(ClientContext &, ReplacementScanInput &input, optional_ptr<ReplacementScanData>)
    -> unique_ptr<TableRef> {
	auto &table_name = input.table_name;
	auto lower_name = StringUtil::Lower(table_name);
	// Check if the table name ends with some common geospatial file extensions
	if (StringUtil::EndsWith(lower_name, ".gpkg") || StringUtil::EndsWith(lower_name, ".fgb")) {

		auto table_function = make_uniq<TableFunctionRef>();
		vector<unique_ptr<ParsedExpression>> children;
		children.push_back(make_uniq<ConstantExpression>(Value(table_name)));
		table_function->function = make_uniq<FunctionExpression>("ST_Read", std::move(children));
		return std::move(table_function);
	}
	// else not something we can replace
	return nullptr;
}

//----------------------------------------------------------------------------------------------------------------------
// REGISTER
//----------------------------------------------------------------------------------------------------------------------
void Register(ExtensionLoader &loader) {
	TableFunction read_func("gdal_read", {LogicalType::VARCHAR}, Scan, Bind, InitGlobal);
	read_func.cardinality = Cardinality;
	read_func.statistics = Statistics;
	read_func.table_scan_progress = Progress;
	read_func.pushdown_complex_filter = Pushdown;

	read_func.named_parameters["open_options"] = LogicalType::LIST(LogicalType::VARCHAR);
	read_func.named_parameters["allowed_drivers"] = LogicalType::LIST(LogicalType::VARCHAR);
	read_func.named_parameters["sibling_files"] = LogicalType::LIST(LogicalType::VARCHAR);
	read_func.named_parameters["layer"] = LogicalType::VARCHAR;
	read_func.named_parameters["max_batch_size"] = LogicalType::INTEGER;
	read_func.named_parameters["keep_wkb"] = LogicalType::BOOLEAN;

	loader.RegisterFunction(read_func);

	auto &config = DBConfig::GetConfig(loader.GetDatabaseInstance());
	config.replacement_scans.emplace_back(ReplacementScan);
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

	CPLStringList driver_options;
	CPLStringList layer_options;

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

bool MatchOption(const char *name, const pair<string, vector<Value>> &option, bool list = false) {
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
				result->layer_options.AddString(val.GetValue<string>().c_str());
			}
			continue;
		}

		if (MatchOption("DATASET_CREATION_OPTIONS", option, true)) {
			for (auto &val : option.second) {
				result->driver_options.AddString(val.GetValue<string>().c_str());
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
		if (srs) {
			OSRDestroySpatialReference(srs);
			srs = nullptr;
		}
	}

	mutex lock;
	GDALDatasetH dataset = nullptr;
	OGRLayerH layer = nullptr;
	OGRSpatialReferenceH srs = nullptr;
};

auto InitGlobal(ClientContext &context, FunctionData &bdata_p, const string &path) -> unique_ptr<GlobalFunctionData> {
	auto &bdata = bdata_p.Cast<BindData>();
	auto result = make_uniq<GlobalState>();

	const auto driver = GDALGetDriverByName(bdata.driver_name.c_str());
	if (!driver) {
		throw InvalidInputException("Could not find GDAL driver: " + bdata.driver_name);
	}

	// Create Dataset
	result->dataset = GDALCreate(driver, bdata.file_path.c_str(), 0, 0, 0, GDT_Unknown, bdata.driver_options);
	if (!result->dataset) {
		throw IOException("Could not create GDAL dataset at: " + bdata.file_path);
	}

	if (!bdata.target_srs.empty()) {
		// Make a new spatial reference object, and set it from the user input
		result->srs = OSRNewSpatialReference(nullptr);
		OSRSetFromUserInput(result->srs, bdata.target_srs.c_str());
	}

	// Create Layer
	result->layer = GDALDatasetCreateLayer(result->dataset, bdata.driver_name.c_str(), result->srs, bdata.geometry_type,
	                                       bdata.layer_options);

	if (!result->layer) {
		throw IOException("Could not create GDAL layer in dataset at: " + bdata.file_path);
	}

	// Create fields for all children
	auto geometry_field_count = 0;
	for (auto i = 0; i < bdata.schema.n_children; i++) {
		const auto child_schema = bdata.schema.children[i];

		// Check if this is a geometry field
		if (child_schema->metadata != nullptr) {
			// TODO: Look for arrow metadata!
			geometry_field_count++;
			if (geometry_field_count > 1) {
				throw NotImplementedException("Multiple geometry fields not supported yet");
			}
		} else {
			// Register normal attribute
			if (!OGR_L_CreateFieldFromArrowSchema(result->layer, child_schema, nullptr)) {
				throw IOException("Could not create field in GDAL layer for column: " + string(child_schema->name));
			}
		}
	}

	return std::move(result);
}

//----------------------------------------------------------------------------------------------------------------------
// Local State
//----------------------------------------------------------------------------------------------------------------------
class LocalState final : public LocalFunctionData {
public:
	~LocalState() override {
		if (array.release) {
			array.release(&array);
			array.release = nullptr;
		}
	}
	ArrowArray array;
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
	auto &lstate = lstate_p.Cast<LocalState>();

	auto &arrow_array = lstate.array;
	auto &arrow_schema = bdata.schema;

	// Convert to Arrow array
	ArrowConverter::ToArrowArray(input, &arrow_array, bdata.props, bdata.extension_type_cast);

	// Sink the Arrow array into GDAL
	{
		// Lock
		lock_guard<mutex> guard(gstate.lock);

		// Sink into GDAL
		OGR_L_WriteArrowBatch(gstate.layer, &arrow_schema, &arrow_array, nullptr);
	}

	// Release the array
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
	// Nothing to do, we don't have any local state that needs to be merged
}

//----------------------------------------------------------------------------------------------------------------------
// Finalize
//----------------------------------------------------------------------------------------------------------------------
void Finalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<GlobalState>();

	// Flush and close the dataset
	GDALFlushCache(gstate.dataset);
	GDALClose(gstate.dataset);
	gstate.dataset = nullptr;
}

CopyFunctionExecutionMode Mode(bool preserve_insertion_order, bool use_batch_index) {
	// Parallel writes have limited utility since we still lock on each write to GDAL layer
	// But in theory we still benefit from the parallel conversion to Arrow arrays, and this also allows
	// the rest of the pipeline to be parallelized if we don't care about insertion order.
	return preserve_insertion_order ? CopyFunctionExecutionMode::REGULAR_COPY_TO_FILE
	                                : CopyFunctionExecutionMode::PARALLEL_COPY_TO_FILE;
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
	info.execution_mode = Mode;
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
