#pragma once

#include <cstddef>

// forward declaration from geos_c.h
struct GEOSGeom_t;
struct GEOSContextHandle_HS;
typedef struct GEOSContextHandle_HS *GEOSContextHandle_t;

namespace duckdb {

class ArenaAllocator;

struct GeosSerde {
	static size_t GetRequiredSize(GEOSContextHandle_t ctx, const GEOSGeom_t *geom);
	static void Serialize(GEOSContextHandle_t ctx, const GEOSGeom_t *geom, char *buffer, size_t buffer_size);
	static GEOSGeom_t *Deserialize(GEOSContextHandle_t ctx, ArenaAllocator &arena, const char *buffer,
	                               size_t buffer_size);
};

} // namespace duckdb
