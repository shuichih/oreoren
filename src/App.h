//
//  App.h
//
//  Created by Shuichi Hayashi on 11/29/11.
//

#ifndef rt_App_h
#define rt_App_h

#include "Config.h"
#include "Timer.h"

class Photon_map;
class BVH;
class IRenderer;
class IView;

class App
{
public:

    App();
    ~App();
    void Run(int argc, const char * argv[]);
    void Render();
    void RenderScene();
    
private:
    bool Init(int argc, const char * argv[]);
    void BuildBVH();
    void ConvertToUint(u8* pColorBuf, Vec3* pRealColorBuf);
    void DrawDebugStuff();
    void DrawBBox();
    void DrawBVH(const IShape* pShape, int level);
    
private:
    
    Config config;
    Photon_map* pPhotonMap_;
    IRenderer* pRenderer_;
    BVH* pBVH_;
    Timer timer_;
    IView* pView_;
    Vec3* pRealColorBuf_;
    u8* pColorBuf_;
};

#endif
