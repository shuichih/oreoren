#ifndef IRenderer_h
#define IRenderer_h

#include "Common.h"

class Config;
class BVH;
class Scene;
class Image;

/**
 * Renderer Interface
 */
class IRenderer
{
protected:
    IRenderer() {};
    
public:
    virtual void SetConfig(const Config& config) = 0;
    virtual void Run(Image& image, const Scene& scene) = 0;
};

#endif
