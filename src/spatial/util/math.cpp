#include "spatial/util/math.hpp"
#include "duckdb/common/string_util.hpp"

namespace duckdb {

#if SPATIAL_USE_GEOS
// We've got this exposed upstream, we just need to wait for the next release
extern "C" int geos_d2sfixed_buffered_n(double f, uint32_t precision, char *result);

template <class T>
static void FormatDouble(T &buffer, double d, int32_t precision) {
	D_ASSERT(precision >= 0 && precision <= 15);
	char buf[512];
	auto len = geos_d2sfixed_buffered_n(d, 15, buf);
	buffer.insert(buffer.end(), buf, buf + len);
}

void MathUtil::format_coord(double x, double y, vector<char> &buffer, int32_t precision) {
	D_ASSERT(precision >= 0 && precision <= 15);
	FormatDouble(buffer, x, precision);
	buffer.push_back(' ');
	FormatDouble(buffer, y, precision);
}

void MathUtil::format_coord(double d, vector<char> &buffer, int32_t precision) {
	FormatDouble(buffer, d, precision);
}

string MathUtil::format_coord(double d) {
	string result;
	FormatDouble(result, d, 15);
	return result;
}

string MathUtil::format_coord(double x, double y) {
	string result;
	FormatDouble(result, x, 15);
	result.push_back(' ');
	FormatDouble(result, y, 15);
	return result;
}

string MathUtil::format_coord(double x, double y, double zm) {
	string result;
	FormatDouble(result, x, 15);
	result.push_back(' ');
	FormatDouble(result, y, 15);
	result.push_back(' ');
	FormatDouble(result, zm, 15);
	return result;
}

string MathUtil::format_coord(double x, double y, double z, double m) {
	string result;
	FormatDouble(result, x, 15);
	result.push_back(' ');
	FormatDouble(result, y, 15);
	result.push_back(' ');
	FormatDouble(result, z, 15);
	result.push_back(' ');
	FormatDouble(result, m, 15);
	return result;
}

#else

void MathUtil::format_coord(double x, double y, vector<char> &buffer, int32_t precision) {
	D_ASSERT(precision >= 0 && precision <= 15);
	auto fmt_str = StringUtil::Format("%%.%df %%.%df", precision, precision);
	auto str = StringUtil::Format(fmt_str, x, y);
	buffer.insert(buffer.end(), str.c_str(), str.c_str() + str.size());
}

void MathUtil::format_coord(double d, vector<char> &buffer, int32_t precision) {
	D_ASSERT(precision >= 0 && precision <= 15);
	auto fmt_str = StringUtil::Format("%%.%df", precision);
	auto str = StringUtil::Format(fmt_str, d);
	buffer.insert(buffer.end(), str.c_str(), str.c_str() + str.size());
}

string MathUtil::format_coord(double d) {
	return StringUtil::Format("%G", d);
}

string MathUtil::format_coord(double x, double y) {
	return StringUtil::Format("%G %G", x, y);
}

string MathUtil::format_coord(double x, double y, double zm) {
	return StringUtil::Format("%G %G %G", x, y, zm);
}

string MathUtil::format_coord(double x, double y, double z, double m) {
	return StringUtil::Format("%G %G %G %G", x, y, z, m);
}

#endif

} // namespace duckdb
