#include "ONB.h"

#define ONB_EPSILON 0.01f

ONB::ONB()
{
}

ONB::ONB(const Vec3& u, const Vec3& v, const Vec3& w)
{
    u_ = u;
    v_ = v;
    w_ = w;
}

void ONB::InitFromW(const Vec3& w)
{
    Vec3 n(1.0f, 0.0f, 0.0f);
    Vec3 m(0.0f, 1.0f, 0.0f);
    
    w_ = w;
    w_.normalize();
    u_ = w_ % n;
    if (u_.length() < ONB_EPSILON) {
        u_ = w_ % m;
    }
    v_ = w_ % u_;
    
}

