#pragma once

#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "spatial/geometry/bbox.hpp"
#include "spatial/geometry/geometry_properties.hpp"
#include "spatial/util/cursor.hpp"

#include "duckdb/common/type_util.hpp"
#include "spatial/util/math.hpp"

namespace duckdb {

enum class LegacyGeometryType : uint8_t {
	POINT = 0,
	LINESTRING,
	POLYGON,
	MULTIPOINT,
	MULTILINESTRING,
	MULTIPOLYGON,
	GEOMETRYCOLLECTION
};

struct LegacyGeometryTypes {
	static bool IsSinglePart(LegacyGeometryType type) {
		return type == LegacyGeometryType::POINT || type == LegacyGeometryType::LINESTRING;
	}

	static bool IsMultiPart(LegacyGeometryType type) {
		return type == LegacyGeometryType::POLYGON || type == LegacyGeometryType::MULTIPOINT ||
		       type == LegacyGeometryType::MULTILINESTRING || type == LegacyGeometryType::MULTIPOLYGON ||
		       type == LegacyGeometryType::GEOMETRYCOLLECTION;
	}

	static bool IsCollection(LegacyGeometryType type) {
		return type == LegacyGeometryType::MULTIPOINT || type == LegacyGeometryType::MULTILINESTRING ||
		       type == LegacyGeometryType::MULTIPOLYGON || type == LegacyGeometryType::GEOMETRYCOLLECTION;
	}

	static string ToString(LegacyGeometryType type) {
		switch (type) {
		case LegacyGeometryType::POINT:
			return "POINT";
		case LegacyGeometryType::LINESTRING:
			return "LINESTRING";
		case LegacyGeometryType::POLYGON:
			return "POLYGON";
		case LegacyGeometryType::MULTIPOINT:
			return "MULTIPOINT";
		case LegacyGeometryType::MULTILINESTRING:
			return "MULTILINESTRING";
		case LegacyGeometryType::MULTIPOLYGON:
			return "MULTIPOLYGON";
		case LegacyGeometryType::GEOMETRYCOLLECTION:
			return "GEOMETRYCOLLECTION";
		default:
			return StringUtil::Format("UNKNOWN(%d)", static_cast<int>(type));
		}
	}
};

enum class SerializedGeometryType : uint32_t {
	POINT = 0,
	LINESTRING,
	POLYGON,
	MULTIPOINT,
	MULTILINESTRING,
	MULTIPOLYGON,
	GEOMETRYCOLLECTION
};

// A serialized geometry
class geometry_t {
private:
	string_t data;

public:
	geometry_t() = default;
	// NOLINTNEXTLINE
	explicit geometry_t(string_t data) : data(data) {
	}

	// NOLINTNEXTLINE
	operator string_t() const {
		return data;
	}

	LegacyGeometryType GetType() const {
		// return the type
		const auto type = Load<LegacyGeometryType>(const_data_ptr_cast(data.GetPrefix()));
		const auto props = Load<GeometryProperties>(const_data_ptr_cast(data.GetPrefix() + 1));
		props.CheckVersion();
		return type;
	}

	GeometryProperties GetProperties() const {
		const auto props = Load<GeometryProperties>(const_data_ptr_cast(data.GetPrefix() + 1));
		// Check the version
		props.CheckVersion();
		return props;
	}

	bool TryGetCachedBounds(Box2D<float> &bbox) const {
		Cursor cursor(data);

		// Read the header
		auto header_type = cursor.Read<LegacyGeometryType>();
		auto properties = cursor.Read<GeometryProperties>();
		auto hash = cursor.Read<uint16_t>();
		(void)hash;

		// Check the version
		properties.CheckVersion();

		if (properties.HasBBox()) {
			cursor.Skip(4); // skip padding

			// Now set the bounding box
			bbox.min.x = cursor.Read<float>();
			bbox.min.y = cursor.Read<float>();
			bbox.max.x = cursor.Read<float>();
			bbox.max.y = cursor.Read<float>();
			return true;
		}

		if (header_type == LegacyGeometryType::POINT) {
			cursor.Skip(4); // skip padding

			// Read the point
			auto type = cursor.Read<SerializedGeometryType>();
			D_ASSERT(type == SerializedGeometryType::POINT);
			(void)type;

			auto count = cursor.Read<uint32_t>();
			if (count == 0) {
				// If the point is empty, there is no bounding box
				return false;
			}

			const auto x = cursor.Read<double>();
			const auto y = cursor.Read<double>();
			bbox.min.x = MathUtil::DoubleToFloatDown(x);
			bbox.min.y = MathUtil::DoubleToFloatDown(y);
			bbox.max.x = MathUtil::DoubleToFloatUp(x);
			bbox.max.y = MathUtil::DoubleToFloatUp(y);
			return true;
		}
		return false;
	}
};

template <>
inline PhysicalType GetTypeId<geometry_t>() {
	return PhysicalType::VARCHAR;
}

static_assert(sizeof(geometry_t) == sizeof(string_t), "geometry_t should be the same size as string_t");

} // namespace duckdb
