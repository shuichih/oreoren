//
//  App.cpp
//
//  Created by Shuichi Hayashi on 11/29/11.
//

#include <algorithm>
#include <cassert>
#include "App.h"
#include "PhotonMap.h"
#include "PhotonMapRenderer.h"
#include "PhotonMapRenderer2.h"
#include "PostEffect.h"
#include "Scene.h"
#include "BVH.h"
#include "RayTracingRenderer.h"
#include "OpenGLView.h"
#include "BmpFileView.h"
#include "Image.h"

using namespace std;

//
//

App::App()
: pBVH_(NULL)
, pImageF32_(NULL)
, pImage_(NULL)
{
}

App::~App()
{
    delete pBVH_;
    delete pImageF32_;
    delete pImage_;
}

void App::Run(int argc, const char * argv[])
{
    if (!Init(argc, argv)) {
        return;
    }

    Render();
    
    pView_->Present(*pImage_);
}

bool App::Init(int argc, const char * argv[])
{
    // Configロード
    const char* path = (argc >= 2) ? argv[1] : "./SceneFiles/render1h.scene";
    if (!config.Load(path)) {
        return false;
    }
    printf("--------------------------------\n\n");
    
    // View初期化
#ifdef __APPLE__
    pView_ = new OpenGLView(&config);
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
    pImageF32_ = new Image(w, h, Image::RGB_F32);
    pImage_ = new Image(w, h, Image::RGBA_8);

    return true;
}

void App::BuildBVH()
{
    if (config.bvhConf.build) {
        Timer timer("Building BVH");
        config.scene.BuildBVH(config.bvhConf.type);
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
    
    pImageF32_->GammaCorrection();
    pImage_->ConvertFrom(*pImageF32_);
}

void App::RenderScene()
{
    // @todo check if light is exist, camera is exist
    // CheckScene();
    
    pRenderer_->Run(*pImageF32_, config.scene);
    
    // Tone Mapping
    if (config.postEffect.toneMapEnabled) {
        ToneMap toneMap;
        toneMap.SetKeyValue(config.postEffect.toneMapKeyValue);
        toneMap.Apply(*pImageF32_, config.windowWidth, config.windowHeight);
    }
}

