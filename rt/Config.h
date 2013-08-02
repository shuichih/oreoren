#ifndef Config_H
#define Config_H

#include "Common.h"
#include "Scene.h"
#include <cstdio>
#include <vector>
#include <string>

using namespace std;

class Material;

//--------------------------------------------------------------------------------
// レンダラ種類
enum RendererType
{
    RTYPE_SIMPLE_RT,
    RTYPE_PHOTON_MAP,
    RTYPE_PHOTON_MAP2
};

//--------------------------------------------------------------------------------
// Section列挙
enum Section
{
    SEC_GENERAL,
    SEC_BVH,
    SEC_RAYTRACING,
    SEC_PMRENDERER,
    SEC_PHOTONMAP,
    SEC_COARSTICPM,
    SEC_SHADOWPM,
    SEC_MATERIAL,
    SEC_CAMERA,
    SEC_LIGHTSOURCE,
    SEC_SPHERE,
    SEC_TRIANGLE,
    SEC_RECTANGLE,
    SEC_CUBOID,
    SEC_SCENEIMPORT,
    SEC_WATERSURFACE,
    SEC_POSTEFFECT,
    SEC_NUM,
};

//--------------------------------------------------------------------------------
// 値の型列挙
enum ItemValueType
{
    IVT_INT,
    IVT_FLOAT,
    IVT_BOOL,
    IVT_VEC2,
    IVT_VEC3,
    IVT_STR,
    IVT_REFL,
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
    Vec3 flux;
    float radius;
    int nSamples;
    
    LightSourceConfig()
    : flux(10000, 10000, 10000)
    , nSamples(64)
    , radius(1)
    {
        typeStr = "POINT";
    }
};

struct SphereConfig
{
    real radius;
    Vec3 position;
    Material* pMaterial;
};

struct SceneImportConfig
{
    string path;
    Vec3 scale;
    Vec3 translate;
    Vec3 rotate;
    Material* pMaterial;
    bool faceReverse;
    
    SceneImportConfig()
    : scale(1, 1, 1)
    , faceReverse(false)
    {}
};

struct PhotonMapConfig
{
    bool enable;
    u32 nPhotons;
    u32 nMaxStorePhotons;
    u32 nEstimatePhotons;
    float estimateDist;
    float estimateEllipseScale;
    u32 maxPhotonBounce;
    bool enableConeFilter;
    float coneFilterK;
    
    PhotonMapConfig()
    : enable(true)
    , nPhotons(100000)
    , nEstimatePhotons(200)
    , estimateDist(15.f)
    , estimateEllipseScale(0.2f)
    , maxPhotonBounce(5)
    , enableConeFilter(true)
    , coneFilterK(1.1f)
    {}
};

struct RayTracingConfig
{
    u32 nSubPixelsSqrt;
    u32 maxRayBounce;
    bool useTentFilter;
    float distanceToProjPlane;
};

struct BVHConfig
{
    bool build;
    BVHType type;
    bool useSIMD;
    bool draw;
    int drawDepth;
    
    BVHConfig()
    : build(true)
    , type(BVH_BINARY)
    , useSIMD(true)
    , draw(false)
    , drawDepth(7)
    {}
};

struct PhotonMapRendererConfig
{
    bool directLight;
    bool indirectLight;
    bool caustic;
    bool shadowEstimate;
    bool drawShadowEstimate;
    u32 nSubPixelsSqrt;
    u32 nTracePhotonsPerThread;
    u32 maxRayBounce;
    bool finalGethering;
    u32 nFinalGetheringRays;
    u32 nMaxGlossyBounce;
    u32 nGlossyRays;
    bool useTentFilter;
    
    PhotonMapRendererConfig()
    : directLight(true)
    , indirectLight(true)
    , caustic(true)
    , shadowEstimate(true)
    , drawShadowEstimate(false)
    , maxRayBounce(5)
    , nTracePhotonsPerThread(10000)
    , useTentFilter(false)
    , finalGethering(false)
    , nFinalGetheringRays(64)
    , nMaxGlossyBounce(1)
    , nGlossyRays(64)
    {}
};

//--------------------------------------------------------------------------------
// セクションパーサ
class SectionParser
{
public:
    SectionParser(const char* pName, Scene* pScene, ItemDesc* paryItemDesc, u32 nItem);
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
    bool ParseVec2(void* pVal, const string& str);
    bool ParseVec3(void* pVal, const string& str);
    bool ParseString(void* pVal, string& str);
    bool ParseRefl(void* pVal, const string& str);
    bool ParseMaterial(void* pVal, const string& str);
    string IntToString(void* pVal);
    string FloatToString(void* pVal);
    string BoolToString(void* pVal);
    string Vec2ToString(void* pVal);
    string Vec3ToString(void* pVal);
    string ReflToString(void* pVal);
    string MaterialToString(void* pVal);
    
protected:
    const char* pName_;
    Scene* pScene_;
    ItemDesc* paryItemDesc_;
    u32 nItem_;
};

//--------------------------------------------------------------------------------
#define SECTION_PARSER_BEGIN(name) \
class name : public SectionParser \
{ \
public: \
    name(const char* pName, Scene* pScene); \
    virtual ~name(); \
    virtual bool OnLeave(); \
private:

#define SECTION_PARSER_END \
};
    
//--------------------------------------------------------------------------------
SECTION_PARSER_BEGIN(MaterialParser)
    Material material_;
SECTION_PARSER_END

//--------------------------------------------------------------------------------
SECTION_PARSER_BEGIN(SphereParser)
    Sphere sphere_;
SECTION_PARSER_END

//--------------------------------------------------------------------------------
SECTION_PARSER_BEGIN(TriangleParser)
    Triangle triangle_;
SECTION_PARSER_END

//--------------------------------------------------------------------------------
SECTION_PARSER_BEGIN(RectangleParser)
    Triangle triangles_[2];
SECTION_PARSER_END

//--------------------------------------------------------------------------------
SECTION_PARSER_BEGIN(CuboidParser)
    Vec3 scale_;
    Vec3 rotate_;
    Vec3 position_;
    Vec3 repeat_;
    int interval_;
    Vec3 margin_;
    bool randColor_;
    Material* pMaterial_;
SECTION_PARSER_END

//--------------------------------------------------------------------------------
SECTION_PARSER_BEGIN(LightSourceParser)
    LightSourceConfig conf_;
    ItemDesc* pItemDesc_;
SECTION_PARSER_END

//--------------------------------------------------------------------------------
SECTION_PARSER_BEGIN(SceneImportParser)
    SceneImportConfig conf_;
    ItemDesc* pItemDesc_;
SECTION_PARSER_END

//--------------------------------------------------------------------------------
SECTION_PARSER_BEGIN(NoiseSurfaceParser)
    Vec3 center_;
    Vec3 scale_;
    Vec3 rotate_;
    Vec2 division_;
    Material* pMaterial_;
    bool noisyHeight_;
    bool noisyColor_;
    ItemDesc* pItemDesc_;
SECTION_PARSER_END

//--------------------------------------------------------------------------------
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
    bool drawBBox;
    RendererType rendererType;
    BVHConfig bvhConf;
    PhotonMapRendererConfig pmRendererConf;
    PhotonMapConfig photonMapConf;
    PhotonMapConfig causticPmConf;
    PhotonMapConfig shadowPmConf;
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
    bool ParseLine(std::string& line);
    bool ParseSceneSection(const string& line);
    
private:
    FILE* fp_;
    char* pBuf_;
    //SectionDesc* pSectionDescs_;
    ItemDesc* pItemDesc_[SEC_NUM];
    //SectionFuncType pCurrFunc;
    //Section currSection;
    SectionParser* parsers_[SEC_NUM];
    SectionParser* pCurrParser_;
    bool bComment_;
    
};

#endif
