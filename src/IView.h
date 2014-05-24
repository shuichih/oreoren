#ifndef IView_h
#define IView_h

#include "Common.h"
#include "IView.h"

class Image;

/**
 * レンダリング結果を表示するView
 */
class IView
{
public:
    virtual bool Init(i32 width, i32 height) = 0;
    virtual bool Present(const Image& image) = 0;
    
protected:
    IView() {}
    virtual ~IView() {}
};

#endif
