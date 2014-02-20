#ifndef _Random_H_
#define _Random_H_

class Random
{
public:
    Random()
    : x_(123456789)
    , y_(362436069)
    , z_(521288629)
    , w_(88675123)
    {}
    
    Random(u32 w)
    : x_(123456789)
    , y_(362436069)
    , z_(521288629)
    , w_(88675123 * w + w)
    {}
    
    void SetSeedW(u32 w)
    {
        w_ = 88675123 * w + w;
    }
    
    u32 U32()
    {
        u32 t = x_ ^ (x_ << 11);
        x_ = y_; y_ = z_; z_ = w_;
        return w_ = (w_ ^ (w_ >> 19)) ^ (t ^ (t >> 8));
    }
    
    float F32()
    {
        return U32() * (1.f / 4294967296.f);
    }
    
private:
    u32 x_;
    u32 y_;
    u32 z_;
    u32 w_;
};


#endif
