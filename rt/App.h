//
//  App.h
//  rt
//
//  Created by 秀一 林 on 11/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef rt_App_h
#define rt_App_h

#include "Config.h"
#include "Timer.h"

class Photon_map;
class BVH;
class IRenderer;

class App
{
public:

    App();
    ~App();
    void Run(int argc, const char * argv[]);
    void Render();
    void RenderScene(Vec3* pRealColorBuf);
    void DrawToBuffer(u8* pRealColorBuf);
    
private:
    void Init(int argc, const char * argv[]);
    void BuildBVH();
    void ConvertToUint(u8* pColorBuf, Vec3* pRealColorBuf);
    void DrawDebugStuff();
    void DrawBBox();
    void DrawBVH(const Shape* pShape, int level);
    
private:
    
    Config config;
    Photon_map* pPhotonMap_;
    IRenderer* pRenderer_;
    BVH* pBVH_;
    Timer timer_;
};

#endif
