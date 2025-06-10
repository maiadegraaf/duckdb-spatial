#define DUCKDB_EXTENSION_MAIN

#include "spatial/spatial_extension.hpp"

#include "duckdb.hpp"
#include "index/rtree/rtree.hpp"
#include "spatial/index/rtree/rtree_module.hpp"
#include "spatial/modules/gdal/gdal_module.hpp"
#if SPATIAL_USE_GEOS
#include "spatial/modules/geos/geos_module.hpp"
#endif
#include "spatial/modules/mvt/mvt_module.hpp"
#include "operators/spatial_operator_extension.hpp"
#include "spatial/modules/main/spatial_functions.hpp"
#include "spatial/modules/osm/osm_module.hpp"
#include "spatial/modules/proj/proj_module.hpp"
#include "spatial/modules/shapefile/shapefile_module.hpp"
#include "spatial/operators/spatial_operator_extension.hpp"
#include "spatial/operators/spatial_join_optimizer.hpp"
#include "spatial/spatial_geoarrow.hpp"
#include "spatial/spatial_types.hpp"

namespace duckdb {

static void LoadInternal(DatabaseInstance &instance) {

	// Register the types
	GeoTypes::Register(instance);

	RegisterSpatialCastFunctions(instance);
	RegisterSpatialScalarFunctions(instance);
	RegisterSpatialAggregateFunctions(instance);
	RegisterSpatialTableFunctions(instance);
	SpatialJoinOptimizer::Register(instance);
	GeoArrow::Register(instance);

	RegisterProjModule(instance);
	RegisterGDALModule(instance);
#if SPATIAL_USE_GEOS
	RegisterGEOSModule(instance);
#endif
	RegisterOSMModule(instance);
	RegisterShapefileModule(instance);
	RegisterMapboxVectorTileModule(instance);

	RTreeModule::RegisterIndex(instance);
	RTreeModule::RegisterIndexPragmas(instance);
	RTreeModule::RegisterIndexScan(instance);
	RTreeModule::RegisterIndexPlanScan(instance);

	RegisterSpatialOperatorExtension(instance);
	;
}

void SpatialExtension::Load(DuckDB &db) {
	LoadInternal(*db.instance);
}

std::string SpatialExtension::Name() {
	return "spatial";
}

} // namespace duckdb

extern "C" {

DUCKDB_EXTENSION_API void spatial_init(duckdb::DatabaseInstance &db) {
	LoadInternal(db);
}

DUCKDB_EXTENSION_API const char *spatial_version() {
	return duckdb::DuckDB::LibraryVersion();
}
}

#ifndef DUCKDB_EXTENSION_MAIN
#error DUCKDB_EXTENSION_MAIN not defined
#endif
