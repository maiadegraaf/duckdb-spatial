#pragma once

namespace duckdb {

class ExtensionLoader;

struct SpatialJoinOptimizer {
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
