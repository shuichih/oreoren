#include "Common.h"

Vec3 operator*(const float f, const Vec3& v)
{
    return Vec3(v.x * f, v.y * f, v.z * f);
}
