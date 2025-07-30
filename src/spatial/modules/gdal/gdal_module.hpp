#pragma once

namespace duckdb {

class ExtensionLoader;

void RegisterGDALModule(ExtensionLoader &loader);

} // namespace duckdb
