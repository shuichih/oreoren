//
//  App.cpp
//
//  Created by Shuichi Hayashi on 11/29/11.
//

#include <algorithm>
#include <cassert>
#include "App.h"
#include "smallpt_fmt.h"
#include "PhotonMap.h"
#include "PhotonMapRenderer.h"
#include "PhotonMapRenderer2.h"
#include "PostEffect.h"
#include "Scene.h"
#include "BVH.h"
#include "RayTracingRenderer.h"
#include "OpenGLView.h"
#include "BmpFileView.h"

using namespace std;

//

inline real clamp(real x)
{
    return x < 0.f ? 0.f : x > 1.f ? 1.f : x;
}

inline int toInt(real x)
{
#ifdef USE_FLOAT
    return int(powf(clamp(x), 1/2.2f) * 255.f + .5f);
#else
    return int(pow(clamp(x), 1/2.2) * 255 + .5);
#endif
}

//

App::App()
: pBVH_(NULL)
, pRealColorBuf_(NULL)
, pColorBuf_(NULL)
{
}

App::~App()
{
    delete pBVH_;
    delete pRealColorBuf_;
    delete pColorBuf_;
}

void App::Run(int argc, const char * argv[])
{
    if (!Init(argc, argv)) {
        return;
    }

    Render();
    
    pView_->Present(pColorBuf_);
}

bool App::Init(int argc, const char * argv[])
{
    // Configロード
    const char* path = (argc >= 2) ? argv[1] : "~/Dev/rt/rt/SceneFiles/simple.scene";
    if (!config.Load(path)) {
        return false;
    }
    printf("--------------------------------\n\n");
    
    // View初期化
#ifdef __APPLE__
    pView_ = new OpenGLView();
#else
    pView_ = new BmpFileView();
#endif
    int w = config.windowWidth;
    int h = config.windowHeight;
    if (!pView_->Init(w, h)) {
        printf("Failed to initialize view.\n");
    }

    // レンダラ初期化
    switch (config.rendererType) {
    case RTYPE_SIMPLE_RT:   pRenderer_ = new RayTracingRenderer(); break;
    case RTYPE_PHOTON_MAP:  pRenderer_ = new PhotonMapRenderer();  break;
    case RTYPE_PHOTON_MAP2: pRenderer_ = new PhotonMapRenderer2(); break;
    }
    pRenderer_->SetConfig(config);

    // 描画バッファ用意
    pRealColorBuf_ = new Vec3[w * h];
    pColorBuf_ = new u8[w * h * sizeof(char)*4];

    return true;
}

void App::BuildBVH()
{
    if (config.bvhConf.build) {
        Timer timer("Building BVH");
        config.scene.BuildBVH(config.bvhConf.type);
    }
}

void App::ConvertToUint(u8* pColorBuf, Vec3* pRealColorBuf)
{
    // Convert to u8 format
    int w = config.windowWidth;
    int h = config.windowHeight;
    for (int i = 0, j=0; i < (w*h); ++i, j+=4)
    {
        pColorBuf[j+0] = (u8)(toInt(pRealColorBuf[i].x)); // including Gamma Correction
        pColorBuf[j+1] = (u8)(toInt(pRealColorBuf[i].y));
        pColorBuf[j+2] = (u8)(toInt(pRealColorBuf[i].z));
        pColorBuf[j+3] = 255;
    }
}

void App::Render()
{
    {
        Timer timer;
        
        BuildBVH();
       
        RenderScene();
        
        timer.PrintElapsed("Total rendering time: ");
    }
    
    ConvertToUint(pColorBuf_, pRealColorBuf_);
}

void App::RenderScene()
{
    // @todo check if light is exist, camera is exist
    // CheckScene();
    
    pRenderer_->Run(pRealColorBuf_, config.scene);
    
    // Tone Mapping
    if (config.postEffect.toneMapEnabled) {
        ToneMap toneMap;
        toneMap.SetKeyValue(config.postEffect.toneMapKeyValue);
        toneMap.Apply(pRealColorBuf_, config.windowWidth, config.windowHeight);
    }
}

