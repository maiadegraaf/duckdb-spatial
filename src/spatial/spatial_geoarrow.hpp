#pragma once

namespace duckdb {

class ExtensionLoader;

struct GeoArrow {
	static void Register(ExtensionLoader &db);
};

} // namespace duckdb