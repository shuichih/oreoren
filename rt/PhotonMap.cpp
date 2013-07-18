// Original code is in
//   Realistic Image Synthesis using Photon Mapping (Japanese Edition)
//   http://ssl.ohmsha.co.jp/cgi-bin/menu.cgi?ISBN=4-274-07950-3

/*
// 使い方の流れ

// photon tracing
foreach (light source) {
	foreach(photon) {
		store(photon);
	}
	scale_photon_power();
}

// irradiance estimate
float irrad[3];
irradiance_estimate(irrad, ...);
irrad *= 1/π	// lambert BRDF

irradiance_estimate()はBRDFを考慮していない。
拡散反射以外の(BRDFが定数でない)場合は、irradiance_estimateに
BRDFを渡してフォトン毎にBRDFを掛ける必要がある。

本の方ではradiance estimateと言っているが多数のフォトン
つまり多数の方向からのradianceを集めて面の色を求めているので
irradiance esitimateの方が正しいと思われる。

// 処理の流れ
balance()
	balance_segment()
		median_split()
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "PhotonMap.h"
#include "PhotonFilter.h"
#include "Common.h"

/**
 * This is the constructor for the photon map.
 * To create the photon map it is necessary to specify the
 * maximum number of photons that will be stored
 */
Photon_map::Photon_map(const int max_phot)
: pFilter_(NULL)
{
	stored_photons = 0;
	half_stored_photons = 0;
	prev_scale = 1;
	max_photons = max_phot;

	photons = (Photon*)malloc(sizeof(Photon) * (max_photons+1)); // 多分[0]はrootなので+1
	if (photons == NULL) {
		fprintf(stderr, "Out of memory initializing photon map¥n");
		exit(-1);
	}
    memset(photons, 0, sizeof(Photon) * (max_photons+1));

	bbox_min[0] = bbox_min[1] = bbox_min[2] = 1e8f;
	bbox_max[0] = bbox_max[1] = bbox_max[2] = -1e8f;

	//----------------------------------------
	// initialize direction conversion tables
	//----------------------------------------

	for (int i=0; i<256; i++) {
		float angle = i*(1.0f/256.0f)*PI;
		costheta[i] = cosf(angle);
		sintheta[i] = sinf(angle);
		cosphi[i]   = cosf(2.0f * angle);
		sinphi[i]   = sinf(2.0f * angle);
	}
}

Photon_map::~Photon_map()
{
	free(photons);
}

void Photon_map::SetFilter(const PhotonFilter* pFilter)
{
    pFilter_ = pFilter;
}

void Photon_map::SetEstimateEllipseScale(float scale)
{
    estimateEllipseScaleInv_ = 1.f / scale;
}
                           
/**
 * photon_dir returns the direction of a photon
 */
void Photon_map::photon_dir(float* dir, const Photon* p) const
{
	dir[0] = sintheta[p->theta] * cosphi[p->phi];
	dir[1] = sintheta[p->theta] * sinphi[p->phi];
	dir[2] = costheta[p->theta];
}

/**
 * irradiance_estimate computes an irradiance estimate
 * at a given surface position
 */
void Photon_map::irradiance_estimate(
	float irrad[3],				// returned irradiance
	const float pos[3],    		// surface position
	const Vec3& normal,   		// surface normal at pos
	const float max_dist,   	// max distance to look for photons
	const int nphotons) const 	// number of photons to use
{
	irrad[0] = irrad[1] = irrad[2] = 0.0f;
	
	// allocaはスタックに確保するので遅くはない。
	NearestPhotons np;
	np.dist2 = (float*)alloca(sizeof(float) * (nphotons+1));
	np.index = (const Photon**)alloca(sizeof(Photon*) * (nphotons+1));
	
	np.pos[0] = pos[0]; np.pos[1] = pos[1]; np.pos[2] = pos[2];
    np.normal = normal;
	np.max = nphotons;
	np.found = 0;
	np.got_heap = 0;
	np.dist2[0] = max_dist * max_dist;

    // locate the nearest photons
    locate_photons(&np, 1);

    // if less than 8 photons return
    if (np.found < 8)
    {
        return;
    }

    float pdir[3];
    float np0_dist = sqrtf(np.dist2[0]);
    
    // sum irradiance from all photons
    for (int i=1; i<=np.found; i++) {
        const Photon* p = np.index[i];
        // the photon_dir call and following if can be omitted (for speed)
        // if the scene does not have any thin surfaces
        // 面の裏にくっついたフォトンを採らないようにしている
        photon_dir(pdir, p);
        if ((pdir[0]*normal.e[0] + pdir[1]*normal.e[1] + pdir[2]*normal.e[2]) < 0.0f) {
            float w = 1.f;
            if (pFilter_) // ループ内で何度もチェックするのが無駄
            {
                float dist = sqrtf(np.dist2[i]);
                w = pFilter_->Weight(dist, np0_dist);
            }
            irrad[0] += p->power[0] * w; // フィルタなしのときw掛けるの無駄
            irrad[1] += p->power[1] * w;
            irrad[2] += p->power[2] * w;
        }
    }

    // BRDFの1/π掛けてないけど、それであってるっぽい
    float tmp = 1.0f / (PI * np.dist2[0]);  // estimate of density
    if (pFilter_)
        tmp *= pFilter_->Normalizer();
    irrad[0] *= tmp;
    irrad[1] *= tmp;
    irrad[2] *= tmp;
}

/**
 * penumbra_estimate computes a ratio of direct light
 * at a given surface position
 */
float Photon_map::penumbra_estimate(
    int lightNo,
	const float pos[3],    		// surface position
	const Vec3& normal,   		// surface normal at pos
	const float max_dist,   	// max distance to look for photons
	const int nphotons) const 	// number of photons to use
{
	// allocaはスタックに確保するので遅くはない。
	NearestPhotons np;
	np.dist2 = (float*)alloca(sizeof(float) * (nphotons+1));
	np.index = (const Photon**)alloca(sizeof(Photon*) * (nphotons+1));
	
	np.pos[0] = pos[0]; np.pos[1] = pos[1]; np.pos[2] = pos[2];
    np.normal = normal;
	np.max = nphotons;
	np.found = 0;
	np.got_heap = 0;
	np.dist2[0] = max_dist * max_dist;

    // locate the nearest photons
    locate_photons(&np, 1);

    // if less than 8 photons return
//    if (np.found < 8)
//    {
//        return 0.5f;
//    }

    float pdir[3];
    
    // 直接光フォトン数と影フォトン数から直接光が当たっている割合を求める
    int nDirect = 0;
    int nShadow = 0;
    for (int i=1; i<=np.found; i++) {
        const Photon* p = np.index[i];
        short photonLightNo = p->flag >> 3;
        if (photonLightNo == lightNo)
        {
            photon_dir(pdir, p);
            if ((pdir[0]*normal.e[0] + pdir[1]*normal.e[1] + pdir[2]*normal.e[2]) < 0.0f) {
                if ((p->flag & Photon::FLAG_DIRECT) != 0) {
                    nDirect++;
                }
                else if (p->power[0] < 0 || p->power[1] < 0 || p->power[2] < 0) {
                    nShadow++;
                }
            }
        }
    }

    if (nDirect+nShadow == 0)
        return 1.0f;

    return (nDirect) / (float)(nDirect+nShadow);
}


/**
 * locate_photons finds the nearest photons in the
 * photon map given the parameters in np
 */
void Photon_map::locate_photons(
    NearestPhotons* const np,
    const int index) const    // 部分木のroot
{
    const Photon* p = &photons[index];
    float dist1;

	// この条件を満たさない場合はもう子がないということ。バランス2分木なのでそうなる
    if (index < half_stored_photons) {
		// フォトンからKD木分割面までの距離
        short plane = (p->flag & 0x0003);
        dist1 = np->pos[plane] - p->pos[plane];

		if (dist1 > 0.0f) { // if dist1 is positive search right plane
			locate_photons(np, 2*index+1);		// ヒープ(整列2分木)に合わせたindex計算
			if (dist1*dist1 < np->dist2[0])
				locate_photons(np, 2*index);
		} else {            // if dist1 is negative search left first
			locate_photons(np, 2*index);
			if (dist1*dist1 < np->dist2[0])
				locate_photons(np, 2*index+1);
		}
	}

	// compute squared distance between current photon and np->pos
	// 推定する点からフォトンへの2乗距離を計算
    Vec3 d(p->pos[0] - np->pos[0], p->pos[1] - np->pos[1], p->pos[2] - np->pos[2]);
    float dist2 = d.lengthSquared();
    
    // 楕円の範囲内のフォトンを集める
    Vec3 nd = np->normal * d.dot(np->normal); // 法線に平行な成分
    Vec3 pn = nd - d; // 法線に垂直な成分
    nd *= estimateEllipseScaleInv_; // 垂直方向にだけスケール
    d = nd + pn;
	float dist2_ellipse = d.lengthSquared();

	//if (dist2_ellipse < np->dist2[0]) {
	if (dist2_ellipse < np->dist2[0]) {
        // we found a photon :) Insert it in the candidate list

		if (np->found < np->max) { // np->maxは集めるフォトンの数
			// heap is not full; use array
			np->found++;
			np->dist2[np->found] = dist2;
			np->index[np->found] = p;
//			np->dist2_ellipse[np->found] = dist2_ellipse;
		} else {
			int j, parent;

			if (np->got_heap == 0) { // Do we need to build the heap?
				// Build heap
				float dst2;
				const Photon* phot;
				int half_found = np->found >> 1;
				for (int k=half_found; k>=1; k--) {
					parent = k;
					phot = np->index[k];
					dst2 = np->dist2[k];
					while (parent <= half_found) {
						j = parent + parent;
						if (j < np->found && np->dist2[j] < np->dist2[j+1])
							j++;
						if (dst2 >= np->dist2[j])
							break;
						np->dist2[parent] = np->dist2[j];
						np->index[parent] = np->index[j];
						parent = j;
					}
					np->dist2[parent] = dst2;
					np->index[parent] = phot;
				}
				np->got_heap = 1;
			}

			// insert new photon into max heap
			// delete largest element, insert new and reorder the heap

			parent = 1;
			j = 2;
			while (j <= np->found) {
				if (j < np->found && np->dist2[j] < np->dist2[j+1])
					j++;
				if (dist2 > np->dist2[j])
					break;
				np->dist2[parent] = np->dist2[j];
				np->index[parent] = np->index[j];
				parent = j;
				j += j;
			}
			if (dist2 < np->dist2[parent]) {
				np->index[parent] = p;
				np->dist2[parent] = dist2;
			}

			np->dist2[0] = np->dist2[1];
		}
	}
}

/*
 * store puts a photon into the flat array that will form
 * the final kd-tree
 *
 * Call this function to store a photon
 */
void Photon_map::store(
	const float power[3],
	const float pos[3],
	const float dir[3],
    bool directLight,
    int lightNo)
{
	if (stored_photons >= max_photons) {
        printf("[WARNING] over max_photons\n");
		return;
    }
    
    #pragma omp flush(stored_photons)
    #pragma omp atomic
	stored_photons++;
    
	Photon* const node = &photons[stored_photons];

	for (int i=0; i<3; i++) {
		node->pos[i] = pos[i];

		if (node->pos[i] < bbox_min[i])
			bbox_min[i] = node->pos[i];
		if (node->pos[i] > bbox_max[i])
			bbox_max[i] = node->pos[i];

		node->power[i] = power[i];
	}

	int theta = int(acosf(dir[2]) * (256.f/PI));
	if (theta > 255)
		node->theta = 255;
	else
		node->theta = (unsigned char)theta;
	
	int phi = int(atan2f(dir[1], dir[0]) * (256.f/(2.f*PI)));
	if (phi > 255)
		node->phi = 255;
	else if (phi < 0)
		node->phi = (unsigned char)(phi + 256);
	else
		node->phi = (unsigned char)phi;
    
    if (directLight)
        node->flag |= Photon::FLAG_DIRECT;
    else
        node->flag &= ~Photon::FLAG_DIRECT;
    
    node->flag |= (lightNo << 3);
}

/**
 * scale_photon_power is used to scale the power of all
 * photons once they have been emitted from the light
 * source. scale = 1 / (#emitted photons)
 * Call this function after each light source is processed.
 */
void Photon_map::scale_photon_power(const float scale)
{
	for (int i=prev_scale; i<= stored_photons; i++) {
		photons[i].power[0] *= scale;
		photons[i].power[1] *= scale;
		photons[i].power[2] *= scale;
	}
	prev_scale = stored_photons + 1;
}

void Photon_map::balance(void)
{
	if (stored_photons > 1) {
		// allocate two temporary arrays for the balancing procedure
		Photon** pa1 = (Photon**)malloc(sizeof(Photon*) * (stored_photons+1));
		Photon** pa2 = (Photon**)malloc(sizeof(Photon*) * (stored_photons+1));

		for (int i=0; i<=stored_photons; i++)
			pa2[i] = &photons[i];

		// pa2にコピーしたフォトン配列をpa１にbalanceして詰める
		balance_segment(pa1, pa2, 1, 1, stored_photons);
		free(pa2);

		// reorganize balanced kd-tree (make a heap)
		int d, j=1, foo=1;
		Photon foo_photon = photons[j];
		
		// pa1のフォトン配列の順に元のphotonsを並べ替える
		for (int i=1; i<=stored_photons; i++) {
			d = (int)(pa1[j] - photons); // photonsにおけるpa1[j]のフォトンのindex
			pa1[j] = NULL;
			if (d != foo)
				photons[j] = photons[d];
			else {
				photons[j] = foo_photon;

				if (i < stored_photons) {
					for (; foo <= stored_photons; foo++)
						if (pa1[foo] != NULL)
							break;
					foo_photon = photons[foo];
					j = foo;
				}
				continue;
			}
			j = d;
		}
		free(pa1);
	}

	half_stored_photons = stored_photons/2 - 1; // half_stored_photonsは数というかindexなので-1
}


#define swap(ph, a, b) { Photon* ph2 = ph[a]; ph[a] = ph[b]; ph[b] = ph2; } 

// median_split splits the photon array into two separate
// piece around the median with all photons below
// the median in the lower half and all photons above
// than the median in the upper half. The comparison
// criteria is the axis (indicated by the axis parameter)
// (inspired by routine in "Algorithm in C++" by sedgewick)
void Photon_map::median_split(
	Photon** p,
	const int start,    	// start of photon block in array
	const int end,			// end of photon block in array
	const int median,		// desired median number
	const int axis)			// axis to split along
{
	int left = start;
	int right = end;

	while (right > left) {
		const float v = p[right]->pos[axis];
		int i = left - 1;
		int j = right;
		for (;;) {
			while (p[++i]->pos[axis] < v)
				;
			while (p[--j]->pos[axis] > v && j > left)
				;
			if (i >= j)
				break;
			swap(p, i, j);
		}

		swap(p, i, right);
		if (i >= median)
			right = i - 1;
		if (i <= median)
			left = i + 1;
	}
}

// See "Realistic image synthesis using Photon Mapping" chapter 6
// for an explanation of this function
void Photon_map::balance_segment(
	Photon** pbal,		// バランス済み配列(out)
	Photon** porg,		// 元の配列
	const int index,
	const int start,
	const int end)
{
	// 2分探索的にmedianのあたりをつける
	int median = 1;
	while ((4 * median) <= (end - start + 1))
		median += median;
	
	if ((3 * median) <= (end - start + 1)) {
		median += median;
		median += start - 1;
	} else
		median = end - median + 1;
	
	//------------------------------------------------
	// find axis to split along
	//------------------------------------------------

	// bboxの一番長い軸を分割軸にする
	int axis = 2;
	if ((bbox_max[0] - bbox_min[0]) > (bbox_max[1] - bbox_min[1]) &&
		(bbox_max[0] - bbox_min[0]) > (bbox_max[2] - bbox_max[2]))
		axis = 0;
	else if ((bbox_max[1] - bbox_min[1]) > (bbox_max[2] - bbox_min[2]))
		axis = 1;

	//------------------------------------------------
	// partition photon block around the median
	//------------------------------------------------

	median_split(porg, start, end, median, axis);
	pbal[index] = porg[median]; // medianのphotonでindex位置を更新
	pbal[index]->flag |= axis;

	//------------------------------------------------
	// recursively balance the left and right block
	//------------------------------------------------

	if (median > start) {
		// balance left segment
		if (start < median - 1) {
			const float tmp = bbox_max[axis];			// bboxの値をバックアップ
			bbox_max[axis] = pbal[index]->pos[axis];	// bboxを狭める
			balance_segment(pbal, porg, 2*index, start, median-1);	// 再帰的に処理
			bbox_max[axis] = tmp;						// bboxの値を復帰
		} else {
			pbal[2*index] = porg[start];
		}
	}

	if (median < end) {
		// balance right segment
		if (median+1 < end) {
			const float tmp = bbox_min[axis];
			bbox_min[axis] = pbal[index]->pos[axis];
			balance_segment(pbal, porg, 2*index+1, median+1, end);
			bbox_min[axis] = tmp;
		} else {
			pbal[2*index+1] = porg[end];
		}
	}
}


