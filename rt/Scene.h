#ifndef _Scene_H_
#define _Scene_H_

#include <math.h>
#include "Common.h"

namespace {
#ifdef USE_FLOAT
    const real EPSILON = 2e-3;
#else
    const real EPSILON = 1e-4;
#endif
}

enum Refl_t
{
    DIFF, SPEC, REFR
}; // material types, used in radiance()

struct Sphere
{
    real rad;       // radius
    Vec p, e, c;    // position, emission, color
    Refl_t refl;    // reflection type (DIFFuse, SPECular, REFRactive)
    
    Sphere(real rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
        rad(rad_), p(p_), e(e_), c(c_), refl(refl_)
    {
    }

    // returns distance, 0 if nohit
    real intersect(const Ray &r) const
    {
        Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        real t;
        real b = op.dot(r.d);
        real det = b * b - op.dot(op) + rad * rad;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);

        return (t=b-det) > EPSILON ? t : ((t=b+det) > EPSILON ? t : 0);
    }
};

#endif
