#ifndef ONB_h
#define ONB_h

#include "Common.h"

//--------------------------------------------------------------------------------
/**
 * Orthonormal Basis
 */
class ONB
{
public:
    ONB();
    ONB(const Vec3& u, const Vec3& v, const Vec3& w);
    
    void InitFromW(const Vec3& w);
    
    Vec3 u_;
    Vec3 v_;
    Vec3 w_;
};

#endif
