#include "Config.h"
#include "PhotonMapRenderer.h"
#include <cstdio>
#include <string>
#include <cstdlib>
#include <cerrno>
#include "StringUtils.h"
#include "MeshLoader.h"
#include "SimpleArithmetic.h"

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

SectionParser::SectionParser(const char* pName, ItemDesc* paryItemDesc, u32 nItem)
: pName_(pName), paryItemDesc_(paryItemDesc), nItem_(nItem)
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
    
    for (int i=0; i<nItem_; i++) {
        ItemDesc& item = paryItemDesc_[i];
        if (key == item.pName) {
            bool result = true;
            switch (item.type) {
            case IVT_INT:   result = ParseInt  (item.pValAddr, val); break;
            case IVT_FLOAT: result = ParseFloat(item.pValAddr, val); break;
            case IVT_BOOL:  result = ParseBool (item.pValAddr, val); break;
            case IVT_VEC3:  result = ParseVec3 (item.pValAddr, val); break;
            case IVT_STR:   result = ParseString(item.pValAddr, val); break;
            case IVT_MTL:   result = ParseMaterial(item.pValAddr, val); break;
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

bool SectionParser::ParseMaterial(void* pVal, const string& str)
{
    Refl_t val;
    if (str == "DIFF") val = DIFF;
    else if (str == "SPEC") val = SPEC;
    else if (str == "REFR") val = REFR;
    else if (str == "PHONGMETAL") val = PHONGMETAL;
    else {
        return false;
    }
    
    *((Refl_t*)pVal) = val;
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
            case IVT_INT:   valStr = IntToString(item.pValAddr);   break;
            case IVT_FLOAT: valStr = FloatToString(item.pValAddr); break;
            case IVT_BOOL:  valStr = BoolToString(item.pValAddr); break;
            case IVT_VEC3:  valStr = Vec3ToString(item.pValAddr);  break;
            case IVT_STR:   valStr = *((string*)item.pValAddr);  break;
            case IVT_MTL:   valStr = MaterialToString(item.pValAddr);  break;
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
    sprintf(str, "%d", *((int*)pVal));
    return str;
}

string SectionParser::FloatToString(void* pVal)
{
    char str[64];
    sprintf(str, "%f", *((float*)pVal));
    return str;
}

string SectionParser::BoolToString(void* pVal)
{
    char str[64];
    sprintf(str, "%s", (*((bool*)pVal)) ? "true" : "false");
    return str;
}

string SectionParser::Vec3ToString(void* pVal)
{
    char str[128];
    Vec3& vec3 = *((Vec3*)pVal);
    sprintf(str, "%f %f %f", vec3.x, vec3.y, vec3.z);
    return str;
}

string SectionParser::MaterialToString(void* pVal)
{
    Refl_t mtl = *((Refl_t*)pVal);
    if (mtl == DIFF) return "DIFF";
    else if (mtl == SPEC) return "SPEC";
    else if (mtl == REFR) return "REFR";
    else if (mtl == PHONGMETAL) return "PHONGMETAL";
    
    return "";
}

//--------------------------------------------------------------------------------

SphereParser::SphereParser(const char* pName, Scene* pScene)
: SectionParser(pName, NULL, 0)
{
    ItemDesc sphereDesc[] = {
        { "radius", IVT_FLOAT, &sphere_.rad },
        { "position", IVT_VEC3, &sphere_.p },
        { "emission", IVT_VEC3, &sphere_.c },
        { "material", IVT_MTL, &sphere_.refl },
    };
    paryItemDesc_ = CreateItemDesc(sphereDesc, ARRAY_SZ(sphereDesc));
    nItem_ = ARRAY_SZ(sphereDesc);
    
    pScene_ = pScene;
}

SphereParser::~SphereParser()
{
}

bool SphereParser::OnLeave()
{
    pScene_->AddShape(new Sphere(sphere_));
    sphere_ = Sphere();
    return true;
}

//--------------------------------------------------------------------------------

TriangleParser::TriangleParser(const char* pName, Scene* pScene)
: SectionParser(pName, NULL, 0)
{
    ItemDesc itemDescs[] = {
        { "p0", IVT_VEC3, &triangle_.p0 },
        { "p1", IVT_VEC3, &triangle_.p1 },
        { "p2", IVT_VEC3, &triangle_.p2 },
        { "emission", IVT_VEC3, &triangle_.color },
        { "material", IVT_MTL, &triangle_.refl },
    };
    paryItemDesc_ = CreateItemDesc(itemDescs, ARRAY_SZ(itemDescs));
    nItem_ = ARRAY_SZ(itemDescs);
    
    pScene_ = pScene;
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

LightSourceParser::LightSourceParser(const char* pName, Scene* pScene)
: SectionParser(pName, NULL, 0)
{
    ItemDesc litSrcDesc[] = {
        { "position", IVT_VEC3, &litSrc_.position_ },
        { "intensity", IVT_VEC3, &litSrc_.intensity_ }
    };
    paryItemDesc_ = CreateItemDesc(litSrcDesc, ARRAY_SZ(litSrcDesc));
    nItem_ = ARRAY_SZ(litSrcDesc);
    
    pScene_ = pScene;
}

LightSourceParser::~LightSourceParser()
{
}

bool LightSourceParser::OnLeave()
{
    pScene_->AddLightSource(new LightSource(litSrc_));
    litSrc_ = LightSource();
    return true;
}

//--------------------------------------------------------------------------------

SceneImportParser::SceneImportParser(const char* pName, Scene* pScene)
: SectionParser(pName, NULL, 0)
{
    ItemDesc itemDesc[] = {
        { "path", IVT_STR, &conf_.path },
        { "scale", IVT_VEC3, &conf_.scale },
        { "translate", IVT_VEC3, &conf_.translate },
        { "material", IVT_MTL, &conf_.material },
        { "faceReverse", IVT_BOOL, &conf_.faceReverse },
    };
    paryItemDesc_ = CreateItemDesc(itemDesc, ARRAY_SZ(itemDesc));
    nItem_ = ARRAY_SZ(itemDesc);
    
    pScene_ = pScene;
}

SceneImportParser::~SceneImportParser()
{
}

bool SceneImportParser::OnLeave()
{
    ObjLoader loader;
    loader.SetFaceReverse(conf_.faceReverse);
    
    Mesh* pMesh = loader.Load(conf_.path.c_str());
    if (pMesh == NULL)
    {
        printf("Obj file load error: %s\n", conf_.path.c_str());
        return false;
    }
    printf("obj file %s is loaded.\n", conf_.path.c_str());
    
    pMesh->scale(conf_.scale);
    pMesh->translate(conf_.translate);
    pScene_->AddShape(pMesh);
    pMesh->CalcBoundingBox();
    pMesh->material_ = conf_.material;
    
    // Bounding BoxとNormal計算
    pMesh->CalcBoundingBox();
    pMesh->CalcFaceNormals();
    
    return true;
}

//--------------------------------------------------------------------------------

Config::Config()
    : windowWidth(256)
    , windowHeight(256)
    , buildBVH(true)
    , drawBBox(false)
    , drawBVH(false)
    , drawBVHDepth(10)
    , rendererType(RTYPE_SIMPLE_RT)
    , pCurrParser_(NULL)
{
    ItemDesc generalDesc[] = {
        { "windowWidth", IVT_INT, &windowWidth },
        { "windowHeight", IVT_INT, &windowHeight },
        { "buildBVH", IVT_BOOL, &buildBVH },
        { "drawBBox", IVT_BOOL, &drawBBox },
        { "drawBVH", IVT_BOOL, &drawBVH },
        { "drawBVHDepth", IVT_INT, &drawBVHDepth },
        { "rendererType", IVT_INT, &rendererType }
    };
    ItemDesc photonMapDesc[] = {
        { "nPhotons", IVT_INT, &photonMapConf.nPhotons },
        { "nEstimatePhotons", IVT_INT, &photonMapConf.nEstimatePhotons },
        { "estimateDist", IVT_FLOAT, &photonMapConf.estimateDist},
        { "estimateEllipseScale", IVT_FLOAT, &photonMapConf.estimateEllipseScale },
        { "nSubPixelSqrt", IVT_INT, &photonMapConf.nSubPixelsSqrt },
        { "enableConeFilter", IVT_BOOL, &photonMapConf.enableConeFilter },
        { "coneFilterK", IVT_FLOAT, &photonMapConf.coneFilterK },
        { "maxPhotonBounce", IVT_INT, &photonMapConf.maxPhotonBounce },
        { "maxRayBounce", IVT_INT, &photonMapConf.maxRayBounce },
        { "useBVH", IVT_BOOL, &photonMapConf.useBVH },
        { "nTracePhotonsPerThread", IVT_INT, &photonMapConf.nTracePhotonsPerThread },
        { "useTentFilter", IVT_BOOL, &photonMapConf.useTentFilter },
        { "distanceToProjPlane", IVT_FLOAT, &photonMapConf.distanceToProjPlane },
        { "finalGethering", IVT_BOOL, &photonMapConf.finalGethering },
        { "nFinalGetheringRays", IVT_INT, &photonMapConf.nFinalGetheringRays },
        { "nMaxGlossyBounce", IVT_INT, &photonMapConf.nMaxGlossyBounce },
        { "nGlossyRays", IVT_INT, &photonMapConf.nGlossyRays }
    };
    ItemDesc rayTracingDesc[] = {
        { "nSubPixelSqrt", IVT_INT, &rayTracingConf.nSubPixelsSqrt },
        { "maxRayBounce", IVT_INT, &rayTracingConf.maxRayBounce },
        { "useBVH", IVT_BOOL, &rayTracingConf.useBVH },
        { "useTentFilter", IVT_BOOL, &rayTracingConf.useTentFilter },
        { "distanceToProjPlane", IVT_FLOAT, &rayTracingConf.distanceToProjPlane }
    };
    ItemDesc postEffectDesc[] = {
        { "toneMap.enabled", IVT_BOOL, &postEffect.toneMapEnabled },
        { "toneMap.keyValue", IVT_FLOAT, &postEffect.toneMapKeyValue }
    };
    ItemDesc cameraDesc[] = {
        { "camera.position", IVT_VEC3, &camera.position },
        { "camera.direction", IVT_VEC3, &camera.direction },
        { "camera.fovY", IVT_FLOAT, &camera.fovY },
    };
    
    ItemDesc* pGeneralItemDesc = CreateItemDesc(generalDesc, ARRAY_SZ(generalDesc));
    ItemDesc* pPhotonMapItemDesc = CreateItemDesc(photonMapDesc, ARRAY_SZ(photonMapDesc));
    ItemDesc* pRayTracingItemDesc = CreateItemDesc(rayTracingDesc, ARRAY_SZ(rayTracingDesc));
    ItemDesc* pPostEffectItemDesc = CreateItemDesc(postEffectDesc, ARRAY_SZ(postEffectDesc));
    ItemDesc* pCameraItemDesc = CreateItemDesc(cameraDesc, ARRAY_SZ(cameraDesc));
 
    parsers_[SEC_GENERAL]     = new SectionParser("[General]",     pGeneralItemDesc, ARRAY_SZ(generalDesc));
    parsers_[SEC_PHOTONMAP]   = new SectionParser("[PhotonMap]",   pPhotonMapItemDesc, ARRAY_SZ(photonMapDesc));
    parsers_[SEC_RAYTRACING]  = new SectionParser("[RayTracing]",  pRayTracingItemDesc, ARRAY_SZ(rayTracingDesc));
    parsers_[SEC_POSTEFFECT]  = new SectionParser("[PostEffect]",  pPostEffectItemDesc, ARRAY_SZ(postEffectDesc));
    parsers_[SEC_CAMERA]      = new SectionParser("[Camera]",      pCameraItemDesc, ARRAY_SZ(cameraDesc));
    parsers_[SEC_LIGHTSOURCE] = new LightSourceParser("[LightSource]", &scene);
    parsers_[SEC_SPHERE]      = new SphereParser("[Sphere]", &scene);
    parsers_[SEC_TRIANGLE]    = new TriangleParser("[Triangle]", &scene);
    parsers_[SEC_SCENEIMPORT] = new SceneImportParser("[SceneImport]", &scene);
}

Config::~Config()
{
    for (int i=0; i<SEC_NUM; i++) {
        delete parsers_[i];
    }
}

bool Config::Load(const char* pPath)
{
    FILE* fp = fopen(pPath, "r");
    if (fp == NULL) {
        printf("%s\n", strerror(errno));
        printf("cannot open the file. %s\n", pPath);
        return false;
    }
    
    const int BUF_SZ = 512;
    char buf[BUF_SZ];
    while (feof(fp) == 0) {
        fgets(buf, BUF_SZ, fp);
        bool stop = !ParseLine(buf);
        if (stop || error_)
            break;
    }
    
    if (pCurrParser_) {
        pCurrParser_->Print();
        pCurrParser_->OnLeave();
    }
    
    fclose(fp);
    return true;
}

bool Config::ParseLine(char* pBuf)
{
    string line = pBuf;
    
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
        printf("[Error] Config: out of section. line=%s\n", pBuf);
        return false;
    }
    
    // セクションの1行をパース
    return pCurrParser_->OnParseLine(line.c_str());
}

