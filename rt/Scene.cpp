#include "Scene.h"
#include <cmath>
#include <cstring>
#include "BBox.h"

//--------------------------------------------------------------------------------

Shape::~Shape()
{
}

//--------------------------------------------------------------------------------

Sphere::Sphere()
: rad(1.f), p(), c(1.f, 1.f, 1.f), refl(DIFF)
{
}

Sphere::Sphere(const Sphere& sphere)
: rad(sphere.rad), p(sphere.p), c(sphere.c), refl(sphere.refl)
{
}

Sphere::Sphere(real rad_, Vec3 p_, Vec3 c_, Refl_t refl_)
 : rad(rad_), p(p_), c(c_), refl(refl_)
{
}

Sphere::~Sphere()
{
}

BBox Sphere::BoundingBox() const
{
    return BBox(p - rad, p + rad);
}

// returns distance, 0 if nohit
bool Sphere::Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const
{
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    Vec3 op = p-r.o;
    real b = op.dot(r.d);
    real det = b * b - op.dot(op) + rad * rad;
    if (det <= 0.f) {
        return false;
    }
    
    det = sqrtf(det);
    
    float t = b - det;
    if (t < tmin || t > tmax) {
        t = b + det;
        if (t < tmin || t > tmax) {
            return false;
        }
    }
    
    rec.t = t;
    rec.normal = ((r.o + r.d * t) - p).normalize();
    rec.color = c;
    rec.refl = refl;
    return true;
}


//--------------------------------------------------------------------------------
Triangle::Triangle(const Vec3& _p0, const Vec3& _p1, const Vec3& _p2, const RGB& _color, Refl_t _refl)
    : p0(_p0), p1(_p1), p2(_p2), color(_color), refl(_refl)
{
    CalcNormal();
}

Triangle::Triangle()
{
    CalcNormal();
}

Triangle::~Triangle()
{
}

BBox Triangle::BoundingBox() const
{
    return BBox(
        std::min(std::min(p0.x, p1.x), p2.x),
        std::min(std::min(p0.y, p1.y), p2.y),
        std::min(std::min(p0.z, p1.z), p2.z),
        std::max(std::max(p0.x, p1.x), p2.x),
        std::max(std::max(p0.y, p1.y), p2.y),
        std::max(std::max(p0.z, p1.z), p2.z)
    );
}

void Triangle::CalcNormal()
{
    normal = ((p1 - p0) % (p2 - p0)).normalize();
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
    Vec3 x = r.o + r.d * t;
    
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
    Vec3 ab = p1 - p0;
    Vec3 ac = p2 - p0;
    Vec3 qp = r.o - (r.o + r.d*100000);
    
    Vec3 n = ab % ac; // 正規化しない
    
    // d < 0なら三角形から離れていく, d == 0なら三角形と平行
    real d = qp.dot(n);
    if (d <= 0.0) return 0;
    
    Vec3 ap = r.o - p0;
    real t = ap.dot(n);
    if (t < 0.0) return false;
    
    // 重心座標の成分を計算し範囲内にあるか判定
    Vec3 e = qp % ap;
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
    rec.normal = Vec3(0, 0, 1);
    rec.color = color;
    rec.refl = refl;
    return true;
}
#endif

#if 1
// based PHISICALLY BASED RENDERING 2ND EDITION, 3.6.2
bool Triangle::Intersect(const Ray& r, float tmin, float tmax, HitRecord& rec) const
{
    Vec3 e1 = p1 - p0;
    Vec3 e2 = p2 - p0;
    Vec3 s1 = r.d % e2;
    real divisor = s1.dot(e1);
    
    if (divisor == 0.0f)
        return false;
    
    real invDivisor = 1.0f / divisor;
    
    // 重心座標を算出して範囲内にあるかチェック
    Vec3 d = r.o - p0;
    real b1 = d.dot(s1) * invDivisor;
    if (b1 < 0.0f || b1 > 1.0f)
        return false;
    
    Vec3 s2 = d % e1;
    real b2 = r.d.dot(s2) * invDivisor;
    if (b2 < 0.0f || b1 + b2 > 1.0f)
        return false;
    
    // t算出
    real t = e2.dot(s2) * invDivisor;
    
    // 光線始点より後ろならヒットしない
    if (t < tmin || t > tmax)
        return false;
    
    rec.t = t;
    rec.normal = normal;
    rec.color = color;
    rec.refl = refl;
    return true;
}
#endif

//--------------------------------------------------------------------------------
MeshTriangle::MeshTriangle() 
    : pMesh(NULL)
{
    indices[0] = 0; 
    indices[1] = 0; 
    indices[2] = 0; 
}

MeshTriangle::~MeshTriangle()
{
}

void MeshTriangle::CalcFaceNormal()
{
    Vec3& p0 = pMesh->pVertices[indices[0]].pos;
    Vec3& p1 = pMesh->pVertices[indices[1]].pos;
    Vec3& p2 = pMesh->pVertices[indices[2]].pos;
    normal = ((p1 - p0) % (p2 - p0)).normalize();
}

BBox MeshTriangle::BoundingBox() const
{
    Vec3 p0 = pMesh->pVertices[indices[0]].pos;
    Vec3 p1 = pMesh->pVertices[indices[1]].pos;
    Vec3 p2 = pMesh->pVertices[indices[2]].pos;
    return BBox(
        std::min(std::min(p0.x, p1.x), p2.x),
        std::min(std::min(p0.y, p1.y), p2.y),
        std::min(std::min(p0.z, p1.z), p2.z),
        std::max(std::max(p0.x, p1.x), p2.x),
        std::max(std::max(p0.y, p1.y), p2.y),
        std::max(std::max(p0.z, p1.z), p2.z)
    );
}

bool MeshTriangle::Intersect(const Ray &r, float tmin, float tmax, HitRecord &rec) const
{
    Vec3 p0 = pMesh->pVertices[indices[0]].pos;
    Vec3 p1 = pMesh->pVertices[indices[1]].pos;
    Vec3 p2 = pMesh->pVertices[indices[2]].pos;
    Vec3 e1 = p1 - p0;
    Vec3 e2 = p2 - p0;
    Vec3 s1 = r.d % e2;
    real divisor = s1.dot(e1);
    
    if (divisor == 0.0f)
        return false;
    
    real invDivisor = 1.0f / divisor;
    
    // 重心座標を算出して範囲内にあるかチェック
    Vec3 d = r.o - p0;
    real b1 = d.dot(s1) * invDivisor;
    if (b1 < 0.0f || b1 > 1.0f)
        return false;
    
    Vec3 s2 = d % e1;
    real b2 = r.d.dot(s2) * invDivisor;
    if (b2 < 0.0f || b1 + b2 > 1.0f)
        return false;
    
    // t算出
    real t = e2.dot(s2) * invDivisor;
    
    // 光線始点より後ろならヒットしない、EPSILONは自己ヒット抑止のため
    if (t < tmin || t > tmax)
        return false;
    
    rec.t = t;
    real b0 = 1.f - b1 - b2;
    rec.normal = ((pMesh->pVertices[indices[0]].normal * b0)
               +  (pMesh->pVertices[indices[1]].normal * b1)
               +  (pMesh->pVertices[indices[2]].normal * b2)).normalize();
    //rec.normal*=-1;
    //rec.normal = normal; // face normal
    
    rec.color = RGB(1, 1, 1);
    rec.refl = pMesh->material_;
    return true;
}

//void MeshTriangle::Reverse()
//{
//    u32 tmp = indices[1];
//    indices[1] = indices[2];
//    indices[2] = tmp;
//}

//--------------------------------------------------------------------------------
Mesh::Mesh(u32 nVertices_, u32 nFaces_)
    : material_(DIFF)
{
    pVertices = new Vertex[nVertices_];
    pFaces = new MeshTriangle[nFaces_];
    ppFaces = new const Shape*[nFaces_];
    nVertices = nVertices_;
    nFaces = nFaces_;
    
    for (int i=0; i<nFaces; i++) {
        pFaces[i].pMesh = this;
        ppFaces[i] = &pFaces[i];
    }
}

Mesh::~Mesh()
{
    delete [] pVertices;
    delete [] pFaces;
}

void Mesh::CalcFaceNormals()
{
    for (int i=0; i<nFaces; i++) {
        pFaces[i].CalcFaceNormal();
    }
}

void Mesh::CalcBoundingBox()
{
    Vec3 minV(FLT_MAX, FLT_MAX, FLT_MAX);
    Vec3 maxV(FLT_MIN, FLT_MIN, FLT_MIN);
    for (u32 i=0; i<nFaces; i++) {
        for (u32 j=0; j<3; j++) {
            u32 index = pFaces[i].indices[j];
            Vec3 p = pVertices[index].pos;
            minV.x = std::min(p.x, minV.x);
            minV.y = std::min(p.y, minV.y);
            minV.z = std::min(p.z, minV.z);
            maxV.x = std::max(p.x, maxV.x);
            maxV.y = std::max(p.y, maxV.y);
            maxV.z = std::max(p.z, maxV.z);
        }
    }
    
    bbox_ = BBox(minV, maxV);
}

BBox Mesh::BoundingBox() const
{
    return bbox_;
}

bool Mesh::Intersect(const Ray &r, float tmin, float tmax, HitRecord& rec) const
{
    bool ret = false;
    for (u32 i=0; i<nFaces; i++) {
        if (pFaces[i].Intersect(r, tmin, tmax, rec)) {
            ret = true;
            tmax = rec.t;
        }
    }
    
    return ret;
}

int Mesh::GetChildNum() const
{
    return nFaces;
    //return 0;
}

const Shape** Mesh::GetChildren() const
{
    return ppFaces;
}

void Mesh::scale(Vec3 scl)
{
    scale(scl.x, scl.y, scl.z);
}

void Mesh::scale(real x, real y, real z)
{
    for (u32 i=0; i<nVertices; i++) {
        pVertices[i].pos.x *= x;
        pVertices[i].pos.y *= y;
        pVertices[i].pos.z *= z;
    }
}

void Mesh::translate(Vec3 transl)
{
    translate(transl.x, transl.y, transl.z);
}

void Mesh::translate(real x, real y, real z)
{
    for (u32 i=0; i<nVertices; i++) {
        pVertices[i].pos.x += x;
        pVertices[i].pos.y += y;
        pVertices[i].pos.z += z;
    }
}

//void Mesh::ReverseFaces()
//{
//    for (u32 i=0; i<nFaces; i++) {
//        pFaces[i].Reverse();
//    }
//}

//--------------------------------------------------------------------------------
Scene::Scene()
{
}

Scene::~Scene()
{
    for (int i=0; i<litSrcs_.size(); i++) {
        delete litSrcs_[i];
    }
    for (int i=0; i<shapes_.size(); i++) {
        delete shapes_[i];
    }
}

void Scene::AddLightSource(const LightSource* pLitSrc)
{
    litSrcs_.push_back(pLitSrc);
}

void Scene::AddShape(const Shape* pShape)
{
    shapes_.push_back(pShape);
}
