#include "Scene.h"
#include <math.h>

namespace {
#ifdef USE_FLOAT
    const real EPSILON = 2e-3;
#else
    const real EPSILON = 2e-4;
#endif
}

//--------------------------------------------------------------------------------

Sphere::Sphere(real rad_, Vec p_, Vec c_, Refl_t refl_)
 : rad(rad_), p(p_), c(c_), refl(refl_)
{
}

// returns distance, 0 if nohit
bool Sphere::intersect(const Ray &r, HitRecord& rec) const
{
    Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    real t;
    real b = op.dot(r.d);
    real det = b * b - op.dot(op) + rad * rad;
    if (det < 0) {
        return false;
    }
    
    det = sqrt(det);
    
    if ((t=b-det) > EPSILON || ((t=b+det) > EPSILON)) {
        rec.t = t;
        rec.normal = ((r.o + r.d * t) - p).norm();
        rec.color = c;
        rec.refl = refl;
        return true;
    }
    
    return false;
}


//--------------------------------------------------------------------------------
Triangle::Triangle(const Vec& _p0, const Vec& _p1, const Vec& _p2, const RGB& _color, Refl_t _refl)
    : p0(_p0), p1(_p1), p2(_p2), color(_color), refl(_refl)
{
    normal = ((p1 - p0) % (p2 - p0)).norm();
}

#if 0
// from Realistic Ray Tracing 2.6.2
// レイを使った交点の式 = 重心座標系での三角形の式 を、クラメルの公式で解く
// でもバグってる。
virtual bool Triangle::intersect(const Ray& r, HitRecord& rec) const
{
    real A = p0.x - p1.x;
    real B = p0.y - p1.y;
    real C = p0.z - p1.z;
    
    real D = p0.x - p2.x;
    real E = p0.y - p2.y;
    real F = p0.z - p2.z;
    
    real G = r.d.x;
    real H = r.d.y;
    real I = r.d.z;
    
    real J = p0.x - r.o.x;
    real K = p0.y - r.o.y;
    real L = p0.z - r.o.z;
    
    real EIHF = E*I-H*F;
    real GFDI = G*F-D*I;
    real DHEG = D*H-E*G;
    
    real denom = (A*EIHF + B*GFDI + C*DHEG);
    real invDenom = 1.0 / denom;

    real beta = (J*EIHF * K*GFDI + L*DHEG) * invDenom;
    
    if (beta <= 0.0 || beta >= 1.0) return false;
    
    real AKJB = A*K - J*B;
    real JCAL = J*C - A*L;
    real BLKC = B*L - K*C;
    
    real gamma = (I*AKJB + H*JCAL + G*BLKC) * invDenom;
    if (gamma <= 0.0 || (beta + gamma) >= 1.0) return false;
    
    real t= -(F*AKJB + E*JCAL + D*BLKC) * invDenom;
    //if (tval >= tmin && tval <= tmax) {// 表と裏のどちらかだけ判定したい場合
    if (t >= EPSILON)
    {
        rec.t = t;
        rec.normal = normal;
        rec.color = color;
        rec.refl = refl;
        return true;
    }
    return false;
}
#endif

#if 0
// 最初にtを算出して、t三角形内にあるか外積で判定する方法。
bool Triangle::intersect(const Ray& r, HitRecord& rec) const
{
    real divisor = r.d.dot(normal);
    if (divisor == 0.0)
        return false;
    
    float t = (p0-r.o).dot(normal) / divisor;
    
    // 光線始点より後ろなら交差しない
    if (t <= EPSILON) return false;
    
    // 光源の向こうなら交差しない（光源の光線方向と交点→ライトへのベクトルが対向しているか見る手もある）
    Vec x = r.o + r.d * t;
    
    // pが各辺の内側にあるか判定
    if (((p1-p0) % (x-p0)).dot(normal) >= 0
     && ((p2-p1) % (x-p1)).dot(normal) >= 0
     && ((p0-p2) % (x-p2)).dot(normal) >= 0)
    {
        rec.t = t;
        rec.normal = normal;
        rec.color = color;
        rec.refl = refl;
        return true;
    }
    return false;
}
#endif

#if 0
// Moller97aの方法。(Real-Time Collision Detection JP Edition, 5.3.6)
// レイの方程式を使った交点の式 = 重心座標系での交点の方程式
// o + td = a + β(b - a) + γ(c - a)
// をクラメルの公式を使って解いている。
// 表しか当たらない。
bool Triangle::intersect(const Ray& r, HitRecord& rec) const
{
    Vec ab = p1 - p0;
    Vec ac = p2 - p0;
    Vec qp = r.o - (r.o + r.d*100000);
    
    Vec n = ab % ac; // 正規化しない
    
    // d < 0なら三角形から離れていく, d == 0なら三角形と平行
    real d = qp.dot(n);
    if (d <= 0.0) return 0;
    
    Vec ap = r.o - p0;
    real t = ap.dot(n);
    if (t < 0.0) return false;
    
    // 重心座標の成分を計算し範囲内にあるか判定
    Vec e = qp % ap;
    real v = ac.dot(e);
    
    if (v < 0.0 || v > d) return false;
    real w = -(ab.dot(e));
    if (w < 0.0 || (v + w) > d) return false;
    
    // tを計算
    real odd = 100000 / d;
    t *= odd;
    
    // 重心座標を計算
    //v *= odd;
    //w *= odd;
    //real u = 1.0 - v - w;
    
    rec.t = t;
    rec.normal = Vec(0, 0, 1);
    rec.color = color;
    rec.refl = refl;
    return true;
}
#endif

#if 1
// based PHISICALLY BASED RENDERING 2ND EDITION, 3.6.2
bool Triangle::intersect(const Ray& r, HitRecord& rec) const
{
    Vec e1 = p1 - p0;
    Vec e2 = p2 - p0;
    Vec s1 = r.d % e2;
    real divisor = s1.dot(e1);
    
    if (divisor == 0.0)
        return false;
    
    real invDivisor = 1.0 / divisor;
    
    // 重心座標を算出して範囲内にあるかチェック
    Vec d = r.o - p0;
    real b1 = d.dot(s1) * invDivisor;
    if (b1 < 0.0 || b1 > 1.0)
        return false;
    
    Vec s2 = d % e1;
    real b2 = r.d.dot(s2) * invDivisor;
    if (b2 < 0.0 || b1 + b2 > 1.0)
        return false;
    
    // t算出
    real t = e2.dot(s2) * invDivisor;
    
    // 光線始点より後ろならヒットしない、EPSILONは自己ヒット抑止のため
    if (t < EPSILON)
        return false;
    
    rec.t = t;
    rec.normal = normal;
    rec.color = color;
    rec.refl = refl;
    return true;
}
#endif

