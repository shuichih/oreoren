#include "Material.h"
#include "Scene.h"
#include "BRDF.h"

//--------------------------------------------------------------------------------
RGB Material::PhotonMapShade(HitRecord& hr)
{
    return RGB(0, 0, 0);
}

//--------------------------------------------------------------------------------
Matte::Matte()
: Material()
, pAmbientBrdf_(new Lambertian())
, pDiffuseBrdf_(new Lambertian())
{
}

void Matte::SetKa(float k)
{
    pAmbientBrdf_->SetKd(k);
}

void Matte::SetKd(float k)
{
    pDiffuseBrdf_->SetKd(k);
}

void Matte::SetCd(const RGB& c)
{
    pAmbientBrdf_->SetCr(c);
    pDiffuseBrdf_->SetCr(c);
}

RGB Matte::PhotonMapShade(HitRecord &hr)
{
#if 0
    pathInfo.diffuseDepth++;
    Vec3 irrad;
    const float BRDF = PI_INV;
    
    // 半影エリアか判定
    float peRatio = .5f; // means penumbra area
    if (pPmRenConf_->shadowEstimate) {
        const float peDist = pShadowPmConf_->estimateDist;
        const int peNum = pShadowPmConf_->nEstimatePhotons;
        peRatio = pShadowPhotonMap_->penumbra_estimate(lightNo_, x.e, nl, peDist, peNum);
    }
    
    // 影度合い別に可視化
    if (pPmRenConf_->drawShadowEstimate) {
        if (peRatio == 1.f) return Vec3(1.f, 1.f, 1.f); // direct
        if (peRatio == 0)   return Vec3(  0,   0,   0); // umbra
        else                return Vec3(.5f, .5f, .5f); // penumbra
    }
    
    // PhotonMappingのLambertではSampleF()でなく
    // F() で帰るBRDFだけ使う
    //pLambert->F(rec, litDir, -r.d);
    //
    // 実際はMaterial::photonmap_shades()でこれを行う。
    // Glossy等さらにトレースする場合もその中で行う
    // diffuseのshadeではShaderRec->world->photonmapと辿れる必要がある
    // PhotonTracingの内容もmaterialに以降する
    // 余裕があればパストレRendererを書く
    
    // HitRecord -> ShaderRec
    // HitRecord.pScene
    // Material::photonmap_shade
    // build new Material in Config
    // Config -> SceneDesc?
    // FIRST: Materialを使ってPMRen::DIFFUSEを処理
    // MaterialをPhotonMapRendererのfriendにするか...とりあえず書いてみよう
    
    // Direct Light
    if (pPmRenConf_->directLight) {
        for (u32 i=0; i<pScene_->GetLightNum(); i++) {
            irrad += BRDF * pScene_->GetLight(i)->DirectLight(x, nl, *pScene_, peRatio, rand);
            assert(irrad.e[0] >= 0);
        }
    }
    
    // Indirect Light
    if (pPmRenConf_->indirectLight) {
        Vec3 tmpIrrad;
        pPhotonMap_->irradiance_estimate(tmpIrrad.e, x.e, nl, pPmConf_->estimateDist, pPmConf_->nEstimatePhotons);
        irrad += BRDF * tmpIrrad;
    }
    
    // Caustics
    if (pPmRenConf_->caustic) {
        Vec3 tmpIrrad;
        pCausticPhotonMap_->irradiance_estimate(tmpIrrad.e, x.e, nl, pCausticPmConf_->estimateDist, pCausticPmConf_->nEstimatePhotons);
        irrad += BRDF * tmpIrrad;
    }
    
    return Vec3(irrad.x * color.x, irrad.y * color.y, irrad.z * color.z);
#endif
    return RGB();
}



