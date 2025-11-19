#include "wkb_module.hpp"

#include "duckdb/common/types.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/common/types/geometry.hpp"

//######################################################################################################################
// Types
//######################################################################################################################
namespace duckdb {
namespace {

struct WKBTypes {

	static LogicalType WKB_BLOB() {
		auto blob_type = LogicalType(LogicalTypeId::BLOB);
		blob_type.SetAlias("WKB_BLOB");
		return blob_type;
	}

	static bool ToWKBCast(Vector &source, Vector &result, idx_t count, CastParameters &) {
		Geometry::ToBinary(source, result, count);
		return true;
	}

	static bool FromWKBCast(Vector &source, Vector &result, idx_t count, CastParameters &params) {
		Geometry::FromBinary(source, result, count, params.strict);
		// TODO: Return false if any errors occurred during the cast
		return true;
	}

	static void Register(ExtensionLoader &loader) {

		// Register the WKB_BLOB type
		loader.RegisterType("WKB_BLOB", WKB_BLOB());

		// Also register casts
		// WKB_BLOB -> GEOMETRY (Explicit)
		loader.RegisterCastFunction(WKB_BLOB(), LogicalType::GEOMETRY(), FromWKBCast);

		// GEOMETRY -> WKB_BLOB (Explicit)
		loader.RegisterCastFunction(LogicalType::GEOMETRY(), WKB_BLOB(), ToWKBCast);

		// WKB_BLOB -> BLOB (Implicit)
		loader.RegisterCastFunction(WKB_BLOB(), LogicalType::BLOB, DefaultCasts::ReinterpretCast, 1);

		// BLOB -> WKB_BLOB (Explicit)
		loader.RegisterCastFunction(LogicalType::BLOB, WKB_BLOB(), DefaultCasts::ReinterpretCast);
	}
};

} // namespace
} // namespace duckdb
//######################################################################################################################
// Functions
//######################################################################################################################
namespace duckdb {
namespace {} // namespace
} // namespace duckdb
//######################################################################################################################
// Module Registration
//######################################################################################################################
namespace duckdb {

void RegisterWKBModule(ExtensionLoader &loader) {
	WKBTypes::Register(loader);
}

} // namespace duckdb
