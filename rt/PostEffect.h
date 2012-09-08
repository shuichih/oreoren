#ifndef PostEffect_h
#define PostEffect_h

#include "Common.h"

/**
 * PostEffect Interface
 */
class IPostEffect
{
public:
    virtual void Apply(Vec3* pBuffer, i32 bufferWidth, i32 bufferHeight) = 0;
};

class ToneMap : public IPostEffect
{
public:
    ToneMap();
    ~ToneMap();
    
    void SetKeyValue(real keyValue);
    void SetDelta(real delta);
    void SetSmallestWhiteLuminance(real l);
    virtual void Apply(Vec3* pBuffer, i32 bufferWidth, i32 bufferHeight);

private:
    real keyValue_;
    real delta_;
    real smallestWhiteLum_;
};

#endif
