#include "Config.h"
#include "PhotonMapRenderer.h"
#include <cstdio>
#include <string>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <cassert>
#include "StringUtils.h"
#include "MeshLoader.h"
#include "SimpleArithmetic.h"
#include "Light.h"
#include "PerlinNoise.h"
#include "OldMaterial.h"
#include "Random.h"
#include "File.h"

using namespace std;

namespace {

ItemDesc* CreateItemDesc(ItemDesc* pInitialDescs, int num)
{
    ItemDesc* pItemDescs = new ItemDesc[num];
    for (int i=0; i<num; i++) {
        pItemDescs[i] = pInitialDescs[i];
    }
    return pItemDescs;
}

};

SectionParser::SectionParser(const char* pName, Scene* pScene, ItemDesc* paryItemDesc, u32 nItem)
: pName_(pName), pScene_(pScene), paryItemDesc_(paryItemDesc), nItem_(nItem)
{
}

SectionParser::~SectionParser()
{
}

bool SectionParser::OnEnter()
{
    return true;
}

bool SectionParser::OnLeave()
{
    return true;
}

bool SectionParser::IsMatch(const string& name)
{
    return name == pName_;
}

bool SectionParser::OnParseLine(const char* pStr)
{
    string line = pStr;
    string key, val;
    if (!ParseKeyValue(key, val, line)) {
        printf("[Error] Config: %s is not Key=Value form.\n", line.c_str());
        return false;
    };
    
    for (u32 i=0; i<nItem_; i++) {
        ItemDesc& item = paryItemDesc_[i];
        if (key == item.pName) {
            bool result = true;
            switch (item.type) {
            case IVT_INT:   result = ParseInt     (item.pValAddr, val); break;
            case IVT_FLOAT: result = ParseFloat   (item.pValAddr, val); break;
            case IVT_BOOL:  result = ParseBool    (item.pValAddr, val); break;
            case IVT_VEC2:  result = ParseVec2    (item.pValAddr, val); break;
            case IVT_VEC3:  result = ParseVec3    (item.pValAddr, val); break;
            case IVT_STR:   result = ParseString  (item.pValAddr, val); break;
            case IVT_REFL:  result = ParseRefl    (item.pValAddr, val); break;
            case IVT_MTL:   result = ParseOldMaterial(item.pValAddr, val); break;
                default:
                result = false;
                break;
            }
            
            if (!result) {
                printf("[Error] Config: parse error. item=%s\n", item.pName);
                return false;
            }
        }
    }
    return true;
}

bool SectionParser::ParseKeyValue(string& key, string& val, const string& line)
{
    vector<string> strs = StringUtils::Split(line, '=');
    if (strs.size() != 2)
        return false;
    
    key = StringUtils::Trim(strs[0]);
    val = StringUtils::Trim(strs[1]);
    return true;
}

bool SectionParser::ParseInt(void* pVal, const string& str)
{
    char* e = NULL;
    long val = strtol(str.c_str(), &e, 10);
    if (e == NULL) {
        return false;
    }
    
    *((int*)pVal) = (int)val;
    
    return true;
}

bool SectionParser::ParseBool(void* pVal, const string& str)
{
    if (StringUtils::Stricmp(str, "true") == 0 || StringUtils::Stricmp(str, "1") == 0) {
        *((bool*)pVal) = true;
    }
    else if (StringUtils::Stricmp(str, "false") == 0 || StringUtils::Stricmp(str, "0") == 0) {
        *((bool*)pVal) = false;
    }
    else {
        return false;
    }
    
    return true;
}


bool SectionParser::ParseFloat(void* pVal, const string& str)
{
    float val;
    if (!SimpleArithmetic::Eval(&val, str)) {
        return false;
    }
    
    *((float*)pVal) = (float)val;
    
    return true;
}

bool SectionParser::ParseVec2(void* pVal, const string& str)
{
    vector<string> strVals = StringUtils::Split(str, ' ');
    if (strVals.size() != 2) {
        return false;
    }
    
    Vec2& vec2 = *((Vec2*)pVal);
    for (int i=0; i < 2; i++) {
        if (!ParseFloat(&vec2.e[i], strVals[i])) {
            return false;
        }
    }
    
    return true;
}

bool SectionParser::ParseVec3(void* pVal, const string& str)
{
    vector<string> strVals = StringUtils::Split(str, ' ');
    if (strVals.size() != 3) {
        return false;
    }
    
    Vec3& vec3 = *((Vec3*)pVal);
    for (int i=0; i < 3; i++) {
        if (!ParseFloat(&vec3.e[i], strVals[i])) {
            return false;
        }
    }
    
    return true;
}

bool SectionParser::ParseString(void* pVal, string& str)
{
    if (str.length() > 0 && str[0] == '\"') {
        str = str.substr(1);
    }
    
    if (str.length() > 0 && str[str.length()-1] == '\"') {
        str = str.substr(0, str.length() - 1);
    }
    *((string*)pVal) = str;
    return true;
}

bool SectionParser::ParseRefl(void* pVal, const string& str)
{
    Refl_t val;
    if (str == "DIFF") val = DIFF;
    else if (str == "SPEC") val = SPEC;
    else if (str == "REFR") val = REFR;
    else if (str == "PHONGMETAL") val = PHONGMETAL;
    else if (str == "LIGHT") val = LIGHT;
    else {
        return false;
    }
    
    *((Refl_t*)pVal) = val;
    return true;
}

bool SectionParser::ParseOldMaterial(void* pVal, const string& str)
{
    OldMaterial* pMtl = pScene_->GetOldMaterial(str);
    if (pMtl) {
        *((OldMaterial**)pVal) = pMtl;
    } else {
        *((OldMaterial**)pVal) = pScene_->GetDefaultOldMaterial();
        printf("WARN!! A material %s is not found.\n", str.c_str());
    }
    return true;
}

void SectionParser::Print()
{
    printf("%s\n", pName_);
    ItemDesc* pItems = paryItemDesc_;
    for (u32 j=0; j<nItem_; j++) {
        ItemDesc& item = pItems[j];
        string valStr;
        switch (item.type) {
            case IVT_INT:   valStr = IntToString(item.pValAddr);      break;
            case IVT_FLOAT: valStr = FloatToString(item.pValAddr);    break;
            case IVT_BOOL:  valStr = BoolToString(item.pValAddr);     break;
            case IVT_VEC2:  valStr = Vec2ToString(item.pValAddr);     break;
            case IVT_VEC3:  valStr = Vec3ToString(item.pValAddr);     break;
            case IVT_STR:   valStr = *((string*)item.pValAddr);       break;
            case IVT_REFL:  valStr = ReflToString(item.pValAddr);     break;
            case IVT_MTL:   valStr = OldMaterialToString(item.pValAddr); break;
            default:
                break;
        }
        printf("%s = %s\n", item.pName, valStr.c_str());
    }
    printf("\n");
}

string SectionParser::IntToString(void* pVal)
{
    char str[64];
    StringUtils::Sprintf(str, sizeof(str), "%d", *((int*)pVal));
    return str;
}

string SectionParser::FloatToString(void* pVal)
{
    char str[64];
    StringUtils::Sprintf(str, sizeof(str), "%f", *((float*)pVal));
    return str;
}

string SectionParser::BoolToString(void* pVal)
{
    char str[64];
    StringUtils::Sprintf(str, sizeof(str), "%s", (*((bool*)pVal)) ? "true" : "false");
    return str;
}

string SectionParser::Vec2ToString(void* pVal)
{
    char str[128];
    Vec2& vec2 = *((Vec2*)pVal);
    StringUtils::Sprintf(str, sizeof(str), "%f %f", vec2.x, vec2.y);
    return str;
}

string SectionParser::Vec3ToString(void* pVal)
{
    char str[128];
    Vec3& vec3 = *((Vec3*)pVal);
    StringUtils::Sprintf(str, sizeof(str), "%f %f %f", vec3.x, vec3.y, vec3.z);
    return str;
}

string SectionParser::ReflToString(void* pVal)
{
    Refl_t refl = *((Refl_t*)pVal);
    if (refl == DIFF) return "DIFF";
    else if (refl == SPEC) return "SPEC";
    else if (refl == REFR) return "REFR";
    else if (refl == PHONGMETAL) return "PHONGMETAL";
    return "";
}

string SectionParser::OldMaterialToString(void* pVal)
{
    OldMaterial* pMtl = *((OldMaterial**)pVal);
    return pMtl->name;
}

//--------------------------------------------------------------------------------
OldMaterialParser::OldMaterialParser(const char* pName, Scene* pScene)
: SectionParser(pName, pScene, NULL, 0)
{
    ItemDesc desc[] = {
        { "name", IVT_STR, &material_.name },
        { "refl", IVT_REFL, &material_.refl },
        { "color", IVT_VEC3, &material_.color },
        { "refrIndex", IVT_FLOAT, &material_.refractiveIndex }
    };
    paryItemDesc_ = CreateItemDesc(desc, ARRAY_SZ(desc));
    nItem_ = ARRAY_SZ(desc);
}

OldMaterialParser::~OldMaterialParser()
{
}

bool OldMaterialParser::OnLeave()
{
    if (pScene_->GetOldMaterial(material_.name)) {
        printf("WARN!! The material %s is already exist.\n", material_.name.c_str());
        return false;
    }
    OldMaterial* pMtl = new OldMaterial(material_);
    pScene_->AddOldMaterial(pMtl->name, pMtl);
    material_.Reset();
    return true;
}

//--------------------------------------------------------------------------------
SphereParser::SphereParser(const char* pName, Scene* pScene)
: SectionParser(pName, pScene, NULL, 0)
, sphere_(1, Vec3(), NULL)
{
    ItemDesc sphereDesc[] = {
        { "radius", IVT_FLOAT, &sphere_.radius },
        { "position", IVT_VEC3, &sphere_.position },
        { "material", IVT_MTL, &sphere_.pOldMaterial },
    };
    paryItemDesc_ = CreateItemDesc(sphereDesc, ARRAY_SZ(sphereDesc));
    nItem_ = ARRAY_SZ(sphereDesc);
}

SphereParser::~SphereParser()
{
}

bool SphereParser::OnLeave()
{
    pScene_->AddShape(new Sphere(sphere_));
    sphere_ = Sphere(1, Vec3(), NULL);
    return true;
}

//--------------------------------------------------------------------------------

TriangleParser::TriangleParser(const char* pName, Scene* pScene)
: SectionParser(pName, pScene, NULL, 0)
{
    ItemDesc itemDescs[] = {
        { "p0", IVT_VEC3, &triangle_.p0 },
        { "p1", IVT_VEC3, &triangle_.p1 },
        { "p2", IVT_VEC3, &triangle_.p2 },
        { "material", IVT_MTL, &triangle_.pOldMaterial },
    };
    paryItemDesc_ = CreateItemDesc(itemDescs, ARRAY_SZ(itemDescs));
    nItem_ = ARRAY_SZ(itemDescs);
}

TriangleParser::~TriangleParser()
{
}

bool TriangleParser::OnLeave()
{
    triangle_.CalcNormal();
    pScene_->AddShape(new Triangle(triangle_));
    triangle_ = Triangle();
    return true;
}


//--------------------------------------------------------------------------------

RectangleParser::RectangleParser(const char* pName, Scene* pScene)
: SectionParser(pName, pScene, NULL, 0)
{
    ItemDesc itemDescs[] = {
        { "p0", IVT_VEC3, &triangles_[0].p0 },
        { "p1", IVT_VEC3, &triangles_[0].p1 },
        { "p2", IVT_VEC3, &triangles_[0].p2 },
        { "p3", IVT_VEC3, &triangles_[1].p2 },
        { "material", IVT_MTL, &triangles_[0].pOldMaterial },
    };
    paryItemDesc_ = CreateItemDesc(itemDescs, ARRAY_SZ(itemDescs));
    nItem_ = ARRAY_SZ(itemDescs);
}

RectangleParser::~RectangleParser()
{
}

bool RectangleParser::OnLeave()
{
    triangles_[1].p0 = triangles_[0].p0;
    triangles_[1].p1 = triangles_[0].p2;
    triangles_[1].pOldMaterial = triangles_[0].pOldMaterial;
    triangles_[0].CalcNormal();
    triangles_[1].CalcNormal();
    pScene_->AddShape(new Triangle(triangles_[0]));
    pScene_->AddShape(new Triangle(triangles_[1]));
    triangles_[0] = Triangle();
    triangles_[1] = Triangle();
    return true;
}

//--------------------------------------------------------------------------------

CuboidParser::CuboidParser(const char* pName, Scene* pScene)
: SectionParser(pName, pScene, NULL, 0)
, scale_(1, 1, 1)
, rotate_()
, position_()
, pOldMaterial_(pScene->GetDefaultOldMaterial())
, repeat_(1, 1, 1)
, interval_(0)
, randColor_(false)
{
    ItemDesc itemDescs[] = {
        { "scale", IVT_VEC3, &scale_ },
        { "rotate", IVT_VEC3, &rotate_ },
        { "position", IVT_VEC3, &position_ },
        { "material", IVT_MTL, &pOldMaterial_ },
        { "repeat", IVT_VEC3, &repeat_ },
        { "interval", IVT_INT, &interval_ },
        { "margin", IVT_VEC3, &margin_ },
        { "randColor", IVT_BOOL, &randColor_ },
    };
    paryItemDesc_ = CreateItemDesc(itemDescs, ARRAY_SZ(itemDescs));
    nItem_ = ARRAY_SZ(itemDescs);
}

CuboidParser::~CuboidParser()
{
}

bool CuboidParser::OnLeave()
{
    // Cubeの頂点と頂点インデックス
    Vec3 v[8] = {
        Vec3(-1, 1,-1),
        Vec3(-1, 1, 1),
        Vec3( 1, 1, 1),
        Vec3( 1, 1,-1),
        Vec3(-1,-1,-1),
        Vec3(-1,-1, 1),
        Vec3( 1,-1, 1),
        Vec3( 1,-1,-1),
    };
    int indices[12][3] = {
        { 0, 1, 2 }, // top
        { 0, 2, 3 },
        { 6, 5, 4 }, // bottom
        { 6, 4, 7 },
        { 0, 4, 5 }, // left
        { 0, 5, 1 },
        { 2, 6, 7 }, // right
        { 2, 7, 3 },
        { 1, 5, 6 }, // front
        { 1, 6, 2 },
        { 3, 7, 4 }, // back
        { 3, 4, 0 },
    };
    
    pScene_->AddOldMaterial(pOldMaterial_->name, pOldMaterial_);
    Random rand;
    rand.SetSeedW(12345678);
    
    // Repeatされた数だけCuboid Meshを作る
    vector<Mesh*> meshes;
    int nRx = int(repeat_.x + 0.5f); // round
    int nRy = int(repeat_.y + 0.5f);
    int nRz = int(repeat_.z + 0.5f);
    int interval = int(interval_ + 0.5f) + 1;
    for (int rx=0; rx<nRx; rx++) {
        for (int ry=0; ry<nRy; ry++) {
            for (int rz=0; rz<nRz; rz++) {
                if ((rx+ry+rz) % interval != 0) continue;
                
                Mesh* pMesh = new Mesh(8, 12, pOldMaterial_);
                Mesh& rMesh = *pMesh;
                meshes.push_back(pMesh);
                
                OldMaterial* pMtl;
                if (randColor_) {
                    float r = rand.F32();
                    float g = rand.F32();
                    float b = rand.F32();
                    const float r1 = 0.3f;
                    const float r2 = 0.5f;
                    r = (r < r1) ? 0.1f : (r > r2) ? 1.0f : 0.6f;
                    g = (g < r1) ? 0.1f : (g > r2) ? 1.0f : 0.6f;
                    b = (b < r1) ? 0.1f : (b > r2) ? 1.0f : 0.6f;
                    //const float r3 = 0.7f;
                    //r = r * (1-r3) + r3;
                    //g = g * (1-r3) + r3;
                    //b = b * (1-r3) + r3;
                    RGB color = Vec3(r, g, b);
                    
                    // 算出された色のマテリアルがなければ作る
                    char mtlName[64];
                    RGB c256 = color * 255.999f;
                    int nRGB [3] = { int(c256.x), int(c256.y), int(c256.z) };
                    StringUtils::Sprintf(mtlName, sizeof(mtlName),
                                         "_D%03d_%03d_%03d_REFL%d_RI%f", nRGB[0], nRGB[1], nRGB[2], pOldMaterial_->refl, pOldMaterial_->refractiveIndex);
                    pMtl = pScene_->GetOldMaterial(mtlName);
                    if (pMtl == NULL) {
                        pMtl = new OldMaterial(mtlName, pOldMaterial_->refl, color, pOldMaterial_->refractiveIndex);
                        pScene_->AddOldMaterial(mtlName, pMtl);
                    }
                } else {
                    pMtl = pOldMaterial_;
                }
                
                for (int i=0; i<8; i++) {
                    rMesh.pVertices[i].pos = v[i];
                }
            
                for (int i=0; i<12; i++) {
                    for (int j=0; j<3; j++) {
                        rMesh.pFaces[i].indices[j] = indices[i][j];
                    }
                    rMesh.pFaces[i].pOldMaterial = pMtl;
                }
                
                float mx = margin_.x;
                float my = margin_.y;
                float mz = margin_.z;
                Vec3 t(
                    -1*(nRx-1 + (nRx-1)*mx) + (rx+rx*mx)*2,
                    -1*(nRy-1 + (nRy-1)*my) + (ry+ry*my)*2,
                    -1*(nRz-1 + (nRz-1)*mz) + (rz+rz*mz)*2
                );
                rMesh.translate(t);
            }
        }
    }
    
    // Meshを統合
    int nMeshes = (int)meshes.size();
    Mesh* pMesh = new Mesh(8*nMeshes, 12*nMeshes, NULL/*@todo*/);
    for (int i=0; i<nMeshes; i++) {
        Mesh& rSrcMesh = *meshes[i];
        int sf = 12*i;
        int sv = 8*i;
        for (u32 j=0; j<rSrcMesh.nVertices; j++) {
            pMesh->pVertices[sv+j] = rSrcMesh.pVertices[j];
        }
        for (u32 j=0; j<rSrcMesh.nFaces; j++) {
            MeshTriangle& rFace = pMesh->pFaces[sf+j];
            rFace = rSrcMesh.pFaces[j];
            rFace.pMesh = pMesh;
            // @todo material
            // rFace.pOldMaterial = pMesh;
            for (int k=0; k<3; k++) {
                pMesh->pFaces[sf+j].indices[k] += sv;
                //printf("%d\n", pMesh->pFaces[s+j].indices[k]);
            }
        }
    }
    for (int i=0; i<nMeshes; i++) {
        delete meshes[i];
    }
    meshes.clear();
    
    // スケール等設定してシーンに追加
    Vec3 scale = scale_ * .5f;
    pMesh->scale(scale);
    pMesh->rotateXYZ(rotate_);
    pMesh->translate(position_);
    pMesh->SetOldMaterial(pOldMaterial_);
    pMesh->CalcBoundingBox();
    pMesh->CalcFaceNormals();
    pMesh->CalcVertexNormals();
    pMesh->SetUseFaceNormal(true);
    pMesh->colorUnit_ = CU_Face;
    pScene_->AddShape(pMesh);
    
    /*
    for (int j=0; j<pMesh->nFaces; j++) {
        for (int k=0; k<3; k++) {
            printf("%d\n", pMesh->pFaces[j].indices[k]);
        }
    }
    */
    
    return true;
}


//--------------------------------------------------------------------------------

LightParser::LightParser(const char* pName, Scene* pScene)
: SectionParser(pName, pScene, NULL, 0)
{
    ItemDesc litDesc[] = {
        // common
        { "type", IVT_STR, &conf_.typeStr },
        { "flux", IVT_VEC3, &conf_.flux },
        // for POINT
        { "position", IVT_VEC3, &conf_.position },
        // for AREA
        { "p0", IVT_VEC3, &conf_.p[0] },
        { "p1", IVT_VEC3, &conf_.p[1] },
        { "p2", IVT_VEC3, &conf_.p[2] },
        { "p3", IVT_VEC3, &conf_.p[3] },
        { "nSamples", IVT_INT, &conf_.nSamples },
        // for SPHERE
        { "radius", IVT_FLOAT, &conf_.radius }
    };
    paryItemDesc_ = CreateItemDesc(litDesc, ARRAY_SZ(litDesc));
    nItem_ = ARRAY_SZ(litDesc);
}

LightParser::~LightParser()
{
}

bool LightParser::OnLeave()
{
    Light* pLitSrc = NULL;
    if (StringUtils::Stricmp(conf_.typeStr, "POINT") == 0) {
        pLitSrc = new PointLight(conf_.position, conf_.flux);
    }
    else if (StringUtils::Stricmp(conf_.typeStr, "AREA") == 0) {
        Vec3* p = conf_.p;
        pLitSrc = new AreaLight(
            p[0], p[1], p[2], p[3],
            conf_.flux,
            conf_.nSamples
        );
        
        OldMaterial* pMtl = pScene_->GetLightOldMaterial();
        Vec3 p2[3] = { p[0], p[2], p[3] };
        IShape* pLitShape0 = new AreaLightShape((AreaLight*)pLitSrc, p, conf_.flux, pMtl);
        IShape* pLitShape1 = new AreaLightShape((AreaLight*)pLitSrc, p2, conf_.flux, pMtl);
        pScene_->AddShape(pLitShape0);
        pScene_->AddShape(pLitShape1);
    }
    else if (StringUtils::Stricmp(conf_.typeStr, "SPHERE") == 0) {
        pLitSrc = new SphereLight(conf_.position, conf_.radius, conf_.flux, conf_.nSamples);
        OldMaterial* pMtl = pScene_->GetLightOldMaterial();
        IShape* pLitShape = new SphereLightShape(
            (SphereLight*)pLitSrc,
            conf_.radius,
            conf_.position,
            conf_.flux,
            pMtl);
        pScene_->AddShape(pLitShape);
    }
    else {
        assert(false);
    }
    pScene_->AddLight(pLitSrc);
                                
    conf_ = LightConfig();
    return true;
}

//--------------------------------------------------------------------------------

SceneImportParser::SceneImportParser(const char* pName, Scene* pScene)
: SectionParser(pName, pScene, NULL, 0)
{
    ItemDesc itemDesc[] = {
        { "path", IVT_STR, &conf_.path },
        { "scale", IVT_VEC3, &conf_.scale },
        { "translate", IVT_VEC3, &conf_.translate },
        { "rotate", IVT_VEC3, &conf_.rotate },
        { "material", IVT_MTL, &conf_.pOldMaterial },
        { "faceReverse", IVT_BOOL, &conf_.faceReverse },
    };
    conf_.pOldMaterial = pScene_->GetDefaultOldMaterial();
    paryItemDesc_ = CreateItemDesc(itemDesc, ARRAY_SZ(itemDesc));
    nItem_ = ARRAY_SZ(itemDesc);
}

SceneImportParser::~SceneImportParser()
{
}

bool SceneImportParser::OnLeave()
{
    MeshLoader loader;
    loader.SetFaceReverse(conf_.faceReverse);
    
    Mesh* pMesh = loader.Load(conf_.path.c_str());
    if (pMesh == NULL)
    {
        printf("Obj file load error: %s\n", conf_.path.c_str());
        return false;
    }
    printf("obj file %s is loaded.\n", conf_.path.c_str());
    
    pMesh->scale(conf_.scale);
    pMesh->rotateXYZ(conf_.rotate);
    pMesh->translate(conf_.translate);
    pScene_->AddShape(pMesh);
    pMesh->CalcBoundingBox();
    pMesh->SetOldMaterial(conf_.pOldMaterial);
    
    // Bounding BoxとNormal計算
    pMesh->CalcBoundingBox();
    //pMesh->CalcFaceNormals();
    
    return true;
}

//--------------------------------------------------------------------------------

NoiseSurfaceParser::NoiseSurfaceParser(const char* pName, Scene* pScene)
: SectionParser(pName, pScene, NULL, 0)
, scale_(1, 1, 1)
, division_(10, 10)
, pOldMaterial_(pScene->GetDefaultOldMaterial())
, noisyHeight_(true)
, noisyColor_(false)
{
    ItemDesc itemDesc[] = {
        { "center", IVT_VEC3, &center_ },
        { "scale", IVT_VEC3, &scale_ },
        { "rotate", IVT_VEC3, &rotate_ },
        { "division", IVT_VEC2, &division_ },
        { "material", IVT_MTL, &pOldMaterial_ },
        { "noisyHeight", IVT_BOOL, &noisyHeight_ },
        { "noisyColor", IVT_BOOL, &noisyColor_ },
    };
    paryItemDesc_ = CreateItemDesc(itemDesc, ARRAY_SZ(itemDesc));
    nItem_ = ARRAY_SZ(itemDesc);
}

NoiseSurfaceParser::~NoiseSurfaceParser()
{
}

bool NoiseSurfaceParser::OnLeave()
{
    PerlinNoise2D noise;
    noise.SetInterporatorType(Interp_Hermite5d);
    
    int divX = (int)division_.e[0];
    int divZ = (int)division_.e[1];
    if (divX <= 0 || divZ <= 0) {
        printf("[Error] Config: divX and divY of a NoiseSurface must be >=1\n");
        return false;
    }
    
    int nVertices = (divX + 1) * (divZ + 1);
    int nFaces = divX * divZ * 2;
    Mesh* pMesh = new Mesh(nVertices, nFaces, pOldMaterial_);
    
    // vertices
    float maxColor = 0;
    float invDivX = (1.f / divX) * 15; // 調整方法考える
    float invDivZ = (1.f / divZ) * 15;
    float amp = 1;
    for (int z=0; z<=divZ; z++) {
        for (int x=0; x<=divX; x++) {
            int iVert = (divX+1) * z + x;
            float y = noisyHeight_
                ? noise.Noise(x*invDivX, z*invDivZ) * amp
                : 0;
            Vec3 pos(
                x/(float)divX - 0.5f,
                y,
                z/(float)divZ - 0.5f
            );
            pMesh->pVertices[iVert].pos = pos;
            if (noisyColor_) {
                Vec3 c = pOldMaterial_->color * noise.Noise(x*invDivX, z*invDivZ) * amp;
                pMesh->pVertices[iVert].color = c;
                maxColor = std::max(maxColor, c.x);
                maxColor = std::max(maxColor, c.y);
                maxColor = std::max(maxColor, c.z);
            }
        }
    }
    if (noisyColor_) {
        for (int i=0; i<nVertices; i++) {
            // normalize vertex color in 0 to 1
            pMesh->pVertices[i].color /= maxColor;
        }
    }
    // faces
    for (int z=0; z<divZ; z++) {
        for (int x=0; x<divX; x++) {
            int iFace = (divX * z + x) * 2;
            int iVerts[4] = {
                (divX+1) * z + x,
                (divX+1) * z + x + 1,
                (divX+1) * (z+1) + x,
                (divX+1) * (z+1) + x + 1,
            };
            pMesh->pFaces[iFace].indices[0] = iVerts[0];
            pMesh->pFaces[iFace].indices[1] = iVerts[2];
            pMesh->pFaces[iFace].indices[2] = iVerts[1];
            pMesh->pFaces[iFace+1].indices[0] = iVerts[1];
            pMesh->pFaces[iFace+1].indices[1] = iVerts[2];
            pMesh->pFaces[iFace+1].indices[2] = iVerts[3];
        }
    }
    
    pMesh->scale(scale_);
    pMesh->rotateXYZ(rotate_);
    pMesh->translate(center_);
    pMesh->SetOldMaterial(pOldMaterial_);
    pMesh->CalcBoundingBox();
    pMesh->CalcFaceNormals();
    pMesh->CalcVertexNormals();
    pMesh->SetUseFaceNormal(false);
    pMesh->colorUnit_ = noisyColor_ ? CU_Vertex : CU_Mesh;
    pMesh->SetOldMaterial(pOldMaterial_);
    pScene_->AddShape(pMesh);
    
    return true;
}

//--------------------------------------------------------------------------------

Config::Config()
    : windowWidth(256)
    , windowHeight(256)
    , drawBBox(false)
    , rendererType(RTYPE_SIMPLE_RT)
    , bgColor(0, 0, 0)
    , pCurrParser_(NULL)
    , bComment_(false)
{
    ItemDesc generalDesc[] = {
        { "windowWidth", IVT_INT, &windowWidth },
        { "windowHeight", IVT_INT, &windowHeight },
        { "drawBBox", IVT_BOOL, &drawBBox },
        { "rendererType", IVT_INT, &rendererType },
        { "bgColor", IVT_VEC3, &bgColor }
    };
    ItemDesc bvhDesc[] = {
        { "build", IVT_BOOL, &bvhConf.build },
        { "type", IVT_INT, &bvhConf.type },
        { "useSIMD", IVT_BOOL, &bvhConf.useSIMD },
        { "draw", IVT_BOOL, &bvhConf.draw },
        { "drawDepth", IVT_INT, &bvhConf.drawDepth },
    };
    ItemDesc pmRendererDesc[] = {
        { "directLight", IVT_BOOL, &pmRendererConf.directLight },
        { "indirectLight", IVT_BOOL, &pmRendererConf.indirectLight },
        { "caustic", IVT_BOOL, &pmRendererConf.caustic },
        { "shadowEstimate", IVT_BOOL, &pmRendererConf.shadowEstimate },
        { "drawShadowEstimate", IVT_BOOL, &pmRendererConf.drawShadowEstimate },
        { "nSubPixelsSqrt", IVT_INT, &pmRendererConf.nSubPixelsSqrt },
        { "maxRayBounce", IVT_INT, &pmRendererConf.maxRayBounce },
        { "nTracePhotonsPerThread", IVT_INT, &pmRendererConf.nTracePhotonsPerThread },
        { "useTentFilter", IVT_BOOL, &pmRendererConf.useTentFilter },
        { "finalGethering", IVT_BOOL, &pmRendererConf.finalGethering },
        { "nFinalGetheringRays", IVT_INT, &pmRendererConf.nFinalGetheringRays },
        { "nMaxGlossyBounce", IVT_INT, &pmRendererConf.nMaxGlossyBounce },
        { "nGlossyRays", IVT_INT, &pmRendererConf.nGlossyRays }
    };
    ItemDesc photonMapDesc[] = {
        { "enable", IVT_BOOL, &photonMapConf.enable },
        { "nPhotons", IVT_INT, &photonMapConf.nPhotons },
        { "nMaxStorePhotons", IVT_INT, &photonMapConf.nMaxStorePhotons },
        { "nEstimatePhotons", IVT_INT, &photonMapConf.nEstimatePhotons },
        { "estimateDist", IVT_FLOAT, &photonMapConf.estimateDist},
        { "estimateEllipseScale", IVT_FLOAT, &photonMapConf.estimateEllipseScale },
        { "enableConeFilter", IVT_BOOL, &photonMapConf.enableConeFilter },
        { "coneFilterK", IVT_FLOAT, &photonMapConf.coneFilterK },
        { "maxPhotonBounce", IVT_INT, &photonMapConf.maxPhotonBounce },
    };
    ItemDesc causticPmDesc[] = {
        { "enable", IVT_BOOL, &causticPmConf.enable },
        { "nPhotons", IVT_INT, &causticPmConf.nPhotons },
        { "nMaxStorePhotons", IVT_INT, &causticPmConf.nMaxStorePhotons },
        { "nEstimatePhotons", IVT_INT, &causticPmConf.nEstimatePhotons },
        { "estimateDist", IVT_FLOAT, &causticPmConf.estimateDist},
        { "estimateEllipseScale", IVT_FLOAT, &causticPmConf.estimateEllipseScale },
        { "enableConeFilter", IVT_BOOL, &causticPmConf.enableConeFilter },
        { "coneFilterK", IVT_FLOAT, &causticPmConf.coneFilterK },
        { "maxPhotonBounce", IVT_INT, &causticPmConf.maxPhotonBounce },
    };
    ItemDesc shadowPmDesc[] = {
        { "enable", IVT_BOOL, &shadowPmConf.enable },
        { "nPhotons", IVT_INT, &shadowPmConf.nPhotons },
        { "nMaxStorePhotons", IVT_INT, &shadowPmConf.nMaxStorePhotons },
        { "nEstimatePhotons", IVT_INT, &shadowPmConf.nEstimatePhotons },
        { "estimateDist", IVT_FLOAT, &shadowPmConf.estimateDist},
        { "estimateEllipseScale", IVT_FLOAT, &shadowPmConf.estimateEllipseScale }
    };
    ItemDesc rayTracingDesc[] = {
        { "nSubPixelsSqrt", IVT_INT, &rayTracingConf.nSubPixelsSqrt },
        { "maxRayBounce", IVT_INT, &rayTracingConf.maxRayBounce },
        { "useTentFilter", IVT_BOOL, &rayTracingConf.useTentFilter },
        { "distanceToProjPlane", IVT_FLOAT, &rayTracingConf.distanceToProjPlane }
    };
    ItemDesc postEffectDesc[] = {
        { "toneMap.enabled", IVT_BOOL, &postEffect.toneMapEnabled },
        { "toneMap.keyValue", IVT_FLOAT, &postEffect.toneMapKeyValue }
    };
    ItemDesc cameraDesc[] = {
        { "position", IVT_VEC3, &camera.position },
        { "direction", IVT_VEC3, &camera.direction },
        { "fovY", IVT_FLOAT, &camera.fovY },
    };
    
    ItemDesc* pGeneralItemDesc = CreateItemDesc(generalDesc, ARRAY_SZ(generalDesc));
    ItemDesc* pBvhItemDesc = CreateItemDesc(bvhDesc, ARRAY_SZ(bvhDesc));
    ItemDesc* pPmRendererItemDesc = CreateItemDesc(pmRendererDesc, ARRAY_SZ(pmRendererDesc));
    ItemDesc* pPhotonMapItemDesc = CreateItemDesc(photonMapDesc, ARRAY_SZ(photonMapDesc));
    ItemDesc* pCausticPmItemDesc = CreateItemDesc(causticPmDesc, ARRAY_SZ(causticPmDesc));
    ItemDesc* pShadowPmItemDesc = CreateItemDesc(shadowPmDesc, ARRAY_SZ(shadowPmDesc));
    ItemDesc* pRayTracingItemDesc = CreateItemDesc(rayTracingDesc, ARRAY_SZ(rayTracingDesc));
    ItemDesc* pPostEffectItemDesc = CreateItemDesc(postEffectDesc, ARRAY_SZ(postEffectDesc));
    ItemDesc* pCameraItemDesc = CreateItemDesc(cameraDesc, ARRAY_SZ(cameraDesc));
 
    parsers_[SEC_GENERAL]     = new SectionParser("[General]", &scene, pGeneralItemDesc, ARRAY_SZ(generalDesc));
    parsers_[SEC_BVH]         = new SectionParser("[BVH]", &scene, pBvhItemDesc, ARRAY_SZ(bvhDesc));
    parsers_[SEC_RAYTRACING]  = new SectionParser("[RayTracing]", &scene, pRayTracingItemDesc, ARRAY_SZ(rayTracingDesc));
    parsers_[SEC_PMRENDERER]  = new SectionParser("[PhotonMapRenderer]", &scene, pPmRendererItemDesc, ARRAY_SZ(pmRendererDesc));
    parsers_[SEC_PHOTONMAP]   = new SectionParser("[PhotonMap]", &scene, pPhotonMapItemDesc, ARRAY_SZ(photonMapDesc));
    parsers_[SEC_COARSTICPM]  = new SectionParser("[CausticPhotonMap]", &scene, pCausticPmItemDesc, ARRAY_SZ(causticPmDesc));
    parsers_[SEC_SHADOWPM]    = new SectionParser("[ShadowPhotonMap]", &scene, pShadowPmItemDesc, ARRAY_SZ(shadowPmDesc));
    parsers_[SEC_MATERIAL]    = new OldMaterialParser("[OldMaterial]", &scene );
    parsers_[SEC_CAMERA]      = new SectionParser("[Camera]", &scene, pCameraItemDesc, ARRAY_SZ(cameraDesc));
    parsers_[SEC_LIGHTSOURCE] = new LightParser("[Light]", &scene);
    parsers_[SEC_SPHERE]      = new SphereParser("[Sphere]", &scene);
    parsers_[SEC_TRIANGLE]    = new TriangleParser("[Triangle]", &scene);
    parsers_[SEC_RECTANGLE]   = new RectangleParser("[Rectangle]", &scene);
    parsers_[SEC_CUBOID]      = new CuboidParser("[Cuboid]", &scene);
    parsers_[SEC_SCENEIMPORT] = new SceneImportParser("[SceneImport]", &scene);
    parsers_[SEC_WATERSURFACE]= new NoiseSurfaceParser("[NoiseSurface]", &scene);
    parsers_[SEC_POSTEFFECT]  = new SectionParser("[PostEffect]", &scene, pPostEffectItemDesc, ARRAY_SZ(postEffectDesc));
}

Config::~Config()
{
    for (int i=0; i<SEC_NUM; i++) {
        delete parsers_[i];
    }
}

bool Config::Load(const char* pPath)
{
    File f;
    if (!f.Open(pPath, "r")) {
        printf("%s\n", f.GetErrorString().c_str());
        printf("cannot open the file. %s\n", pPath);
        return false;
    }
    
    while (!f.IsEof()) {
        std::string line = f.GetLineS();
        bool stop = !ParseLine(line);
        if (stop)
            break;
    }
    
    if (pCurrParser_) {
        pCurrParser_->Print();
        pCurrParser_->OnLeave();
    }
    
    return true;
}

bool Config::ParseLine(std::string& line)
{
    const char* pLine = line.c_str();
    if (pLine[0] == (char)0xEF && pLine[1] == (char)0xBB && pLine[2] == (char)0xBF) {
        line = &pLine[3];
    }
    
    string::size_type acs = line.find("/#");
    if (acs != string::npos) {
        bComment_ = true;
        return true;
    }
    string::size_type ace = line.find("#/");
    if (ace != string::npos) {
        bComment_ = false;
        return true;
    }
    
    if (bComment_)
    {
        return true;
    }
    
    // コメント部分削除
    string::size_type pos = line.find('#');
    if (pos != string::npos) {
        line = line.substr(0, pos);
    }
    
    line = StringUtils::Trim(line);
    
    if (line.length() == 0)
        return true;
    
    // セクション選ぶ
    for (int i=0; i<SEC_NUM; i++) {
        if (parsers_[i]->IsMatch(line)) {
            if (pCurrParser_) {
                pCurrParser_->Print();
                pCurrParser_->OnLeave();
            }
            
            pCurrParser_ = parsers_[i];
            pCurrParser_->OnEnter();
            
            return true;
        }
    }
    if (pCurrParser_ == NULL) {
        printf("[Error] Config: out of section. line=%s\n", line.c_str());
        return false;
    }
    
    // セクションの1行をパース
    return pCurrParser_->OnParseLine(line.c_str());
}

