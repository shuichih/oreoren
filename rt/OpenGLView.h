#ifndef Rayzer_OpenGLView_H
#define Rayzer_OpenGLView_H

#include "Common.h"
#include "IView.h"

#ifdef __APPLE__

class IShape;
class Config;
class BVH;

/**
 * OpenGL Windowに表示するView
 */
class OpenGLView : public IView
{
public:
    
    OpenGLView(Config* pConfig, BVH* pBVH);
    ~OpenGLView();
    virtual bool Init(i32 width, i32 height);
    virtual bool Present(u8* pColorBuf);
    bool DrawFrame();
    
private:
    void DrawToBuffer(u8* pColorBuf);
    void DrawDebugStuff();
    void DrawBVH(const IShape* pShape, int depth);
    void DrawBBox();

    bool initialized_;
    int width_;
    int height_;
    u8* pColorBuf_;
    Config* pConf_;
    BVH* pBVH_;
};

#endif // __APPLE__

#endif // Rayzer_OpenGLView_H
