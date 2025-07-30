#pragma once

namespace duckdb {

class ExtensionLoader;

void RegisterMapboxVectorTileModule(ExtensionLoader &loader);

} // namespace duckdb
