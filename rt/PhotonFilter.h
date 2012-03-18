#ifndef _PhotonFilter_H_
#define _PhotonFilter_H_

// 放射輝度推定する位置とフォトンの距離に応じてフォトンの寄与を調節するフィルタのインターフェース
class PhotonFilter
{
public:
    virtual ~PhotonFilter() {};
    
    virtual float Weight(float dist, float maxDist) const = 0;
    virtual float Normalizer() const = 0;
};

// 円錐フィルタ
class ConeFilter : public PhotonFilter
{
public:
    // kは1.0以上の定数で、1.0に近いほどフィルタ効果が強い
    ConeFilter(float k);
    
    virtual float Weight(float dist, float maxDist) const;
    virtual float Normalizer() const;

private:
    float k_;
    float normalizer_;
};


#endif
