#ifndef Rayzer_OpenGLView_H
#define Rayzer_OpenGLView_H

#include "Common.h"
#include "IView.h"

#ifdef __APPLE__

/**
 * OpenGL Windowに表示するView
 */
class OpenGLView : public IView
{
public:
    
    virtual bool Init(i32 width, i32 height);
    virtual bool Present(u8* pColorBuf);

private:
    bool DrawFrame();
    void DrawToBuffer(u8* pColorBuf);
    void DrawDebugStuff();
    void DrawBVH(const IShape* pShape, int depth);
    void DrawBBox();

    bool initialized_;
    int width_;
    int height_;
};

#endif // __APPLE__

#endif // Rayzer_OpenGLView_H
