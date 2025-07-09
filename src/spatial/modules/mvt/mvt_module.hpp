#pragma once

namespace duckdb {

class DatabaseInstance;

void RegisterMapboxVectorTileModule(DatabaseInstance &db);

} // namespace duckdb