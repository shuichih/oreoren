#include "ONB.h"

#define ONB_EPSILON 0.01f

ONB::ONB()
{
}

ONB::ONB(const Vec3& u__, const Vec3& v__, const Vec3& w__)
{
    u = u__;
    v = v__;
    w = w__;
}

ONB ONB::InitFromW(const Vec3& w__)
{
    Vec3 n(1.0f, 0.0f, 0.0f);
    Vec3 m(0.0f, 1.0f, 0.0f);
    
    Vec3 w = w__;
    w.normalize();
    Vec3 u = w ^ n;
    if (u.length() < ONB_EPSILON) {
        u = w ^ m;
    }
    Vec3 v = w ^ u;
    
    return ONB(u, v, w);
}

