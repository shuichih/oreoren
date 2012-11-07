#ifndef Config_H
#define Config_H

#include "Common.h"
#include "Scene.h"
#include <cstdio>
#include <vector>
#include <string>

using namespace std;

//--------------------------------------------------------------------------------
// レンダラ種類
enum RendererType
{
    RTYPE_SIMPLE_RT,
    RTYPE_PHOTON_MAP
};

//--------------------------------------------------------------------------------
// Section列挙
enum Section
{
    SEC_GENERAL,
    SEC_PHOTONMAP,
    SEC_RAYTRACING,
    SEC_POSTEFFECT,
    SEC_CAMERA,
    SEC_LIGHTSOURCE,
    SEC_SPHERE,
    SEC_TRIANGLE,
    SEC_RECTANGLE,
    SEC_SCENEIMPORT,
    SEC_NUM,
};

//--------------------------------------------------------------------------------
// 値の型列挙
enum ItemValueType
{
    IVT_INT,
    IVT_FLOAT,
    IVT_BOOL,
    IVT_VEC3,
    IVT_STR,
    IVT_MTL,
};

//--------------------------------------------------------------------------------
// 設定アイテム1つ分の定義情報
struct ItemDesc
{
    const char* pName;
    ItemValueType type;
    void* pValAddr;
};


//--------------------------------------------------------------------------------
// 種類別Config
struct CameraConfig
{
    Vec3 position;
    Vec3 direction;
    real fovY;
    
    CameraConfig()
    : direction(0, 0, -1)
    , fovY(cosf(60))
    {}
};
    
struct PostEffectConfig
{
    bool toneMapEnabled;
    real toneMapKeyValue;
    
    PostEffectConfig()
    : toneMapEnabled(false)
    , toneMapKeyValue(0.045f)
    {}
};

struct LightSourceConfig
{
    string typeStr;
    Vec3 p[4];
    Vec3 position;
    Vec3 intensity;
    
    LightSourceConfig()
    : intensity(10000, 10000, 10000)
    {
        typeStr[0] = 'P';
        typeStr[1] = 'O';
        typeStr[2] = 'I';
        typeStr[3] = 'N';
        typeStr[4] = 'T';
    }
};

struct SphereConfig
{
    real radius;
    Vec3 position;
    Vec3 emission;
    Refl_t material;
};

struct SceneImportConfig
{
    string path;
    Vec3 scale;
    Vec3 translate;
    Refl_t material;
    bool faceReverse;
    Vec3 color;
    
    SceneImportConfig()
    : scale(1, 1, 1)
    , faceReverse(false)
    {}
};

struct PhotonMapConfig
{
    u32 nSubPixelsSqrt;
    u32 nPhotons;
    u32 nEstimatePhotons;
    float estimateDist;
    float estimateEllipseScale;
    bool enableConeFilter;
    float coneFilterK;
    u32 maxPhotonBounce;
    u32 maxRayBounce;
    bool useBVH;
    u32 nTracePhotonsPerThread;
    bool useTentFilter;
    float distanceToProjPlane;
    bool finalGethering;
    u32 nFinalGetheringRays;
    u32 nMaxGlossyBounce;
    u32 nGlossyRays;
    
    PhotonMapConfig()
    : nSubPixelsSqrt(1)
    , nPhotons(100000)
    , nEstimatePhotons(200)
    , estimateDist(15.f)
    , estimateEllipseScale(0.2f)
    , enableConeFilter(true)
    , coneFilterK(1.1f)
    , maxPhotonBounce(5)
    , maxRayBounce(5)
    , useBVH(true)
    , nTracePhotonsPerThread(10000)
    , useTentFilter(false)
    , distanceToProjPlane(140)
    , finalGethering(false)
    , nFinalGetheringRays(64)
    , nMaxGlossyBounce(1)
    , nGlossyRays(64)
    {
    }
};

struct RayTracingConfig
{
    u32 nSubPixelsSqrt;
    u32 maxRayBounce;
    bool useBVH;
    bool useTentFilter;
    float distanceToProjPlane;
};


//--------------------------------------------------------------------------------
// セクションパーサ
class SectionParser
{
public:
    SectionParser(const char* pName, ItemDesc* paryItemDesc, u32 nItem);
    virtual ~SectionParser();
    
    virtual bool OnEnter();
    virtual bool OnParseLine(const char* pStr);
    virtual bool OnLeave();
    virtual void Print();
    bool IsMatch(const string& name);
    bool ParseSection(const string& line);
    bool ParseKeyValue(string& key, string& val, const string& line);
    bool ParseInt(void* pVal, const string& val);
    bool ParseFloat(void* pVal, const string& str);
    bool ParseBool(void* pVal, const string& str);
    bool ParseVec3(void* pVal, const string& str);
    bool ParseString(void* pVal, string& str);
    bool ParseMaterial(void* pVal, const string& str);
    string IntToString(void* pVal);
    string FloatToString(void* pVal);
    string BoolToString(void* pVal);
    string Vec3ToString(void* pVal);
    string MaterialToString(void* pVal);
    
    
protected:
    const char* pName_;
    ItemDesc* paryItemDesc_;
    u32 nItem_;
};

//--------------------------------------------------------------------------------
// Sphereセクションパーサ
class SphereParser : public SectionParser
{
public:
    SphereParser(const char* pName, Scene* pScene);
    virtual ~SphereParser();
    
    virtual bool OnLeave();
    
private:
    Sphere sphere_;
    Scene* pScene_;
};

//--------------------------------------------------------------------------------
// Triangle
class TriangleParser : public SectionParser
{
public:
    TriangleParser(const char* pName, Scene* pScene);
    virtual ~TriangleParser();
    
    virtual bool OnLeave();
    
private:
    Triangle triangle_;
    Scene* pScene_;
};

//--------------------------------------------------------------------------------
// Rectangleセクションパーサ
class RectangleParser : public SectionParser
{
public:
    RectangleParser(const char* pName, Scene* pScene);
    virtual ~RectangleParser();
    
    virtual bool OnLeave();
    
private:
    Triangle triangles_[2];
    Scene* pScene_;
};


//--------------------------------------------------------------------------------
// LightSourceセクションパーサ
class LightSourceParser : public SectionParser
{
public:
    LightSourceParser(const char* pName, Scene* pScene);
    virtual ~LightSourceParser();
    
    virtual bool OnLeave();
    
private:
    LightSourceConfig conf_;
    
    ItemDesc* pItemDesc_;
    Scene* pScene_;
};

//--------------------------------------------------------------------------------
// SceneImportセクションパーサ
class SceneImportParser : public SectionParser
{
public:
    SceneImportParser(const char* pName, Scene* pScene);
    virtual ~SceneImportParser();
    
    virtual bool OnLeave();
    
private:
    SceneImportConfig conf_;
    ItemDesc* pItemDesc_;
    Scene* pScene_;
};


//--------------------------------------------------------------------------------
// Config
class Config
{
public:
    Config();
    ~Config();
    bool Load(const char* pPath);
    void Print();
    
public:
    int windowWidth;
    int windowHeight;
    bool buildBVH;
    bool drawBVH;
    bool drawBBox;
    int drawBVHDepth;
    RendererType rendererType;
    PhotonMapConfig photonMapConf;
    RayTracingConfig rayTracingConf;
    PostEffectConfig postEffect;
    CameraConfig camera;
    LightSourceConfig lightSource;
    SphereConfig sphere;
    SceneImportConfig sceneImport;
    Scene scene;
    
    // 以下パース処理用
private:
    
    // セクションを処理するメンバ関数の関数ポインタ型
    typedef bool (Config::*SectionFuncType)(const string& str);
    
private:
    bool ParseLine(char* pBuf);
    bool ParseSceneSection(const string& line);
    
private:
    FILE* fp_;
    char* pBuf_;
    bool error_;
    //SectionDesc* pSectionDescs_;
    ItemDesc* pItemDesc_[SEC_NUM];
    //SectionFuncType pCurrFunc;
    //Section currSection;
    SectionParser* parsers_[SEC_NUM];
    SectionParser* pCurrParser_;
    
};

#endif
