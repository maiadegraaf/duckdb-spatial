#pragma once

namespace duckdb {

class ExtensionLoader;

void RegisterProjModule(ExtensionLoader &db);

} // namespace duckdb