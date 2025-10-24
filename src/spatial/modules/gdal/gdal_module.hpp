#pragma once
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class ExtensionLoader;

void RegisterGDALModule(ExtensionLoader &loader);
void RegisterExtraFunction(ExtensionLoader &loader);
} // namespace duckdb
