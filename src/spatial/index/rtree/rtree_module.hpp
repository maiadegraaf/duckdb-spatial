#pragma once

namespace duckdb {

class ExtensionLoader;

struct RTreeModule {
	static void RegisterIndex(ExtensionLoader &loader);
	static void RegisterIndexScan(ExtensionLoader &loader);
	static void RegisterIndexPlanScan(ExtensionLoader &loader);
	static void RegisterIndexPragmas(ExtensionLoader &loader);
};

} // namespace duckdb
