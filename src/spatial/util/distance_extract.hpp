#pragma once

namespace duckdb {

struct FunctionData;

// Helper class to extract the constant distance from the ST_DWithin function.
struct ST_DWithinHelper {
	static bool TryGetConstDistance(const unique_ptr<FunctionData> &bind_data, double &result);
};

} // namespace duckdb