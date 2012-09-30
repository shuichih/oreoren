#ifndef IRenderer_h
#define IRenderer_h

class Config;
class BVH;
class Scene;
class Vec3;

/**
 * Renderer Interface
 */
class IRenderer
{
protected:
    IRenderer() {};
    
public:
    virtual void SetConfig(const Config& config) = 0;
    virtual void Run(Vec3* pColorBuf, const Scene& scene, BVH* pBVH) = 0;
};

#endif
