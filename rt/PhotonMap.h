// Original code is in
//   Realistic Image Synthesis using Photon Mapping (Japanese Edition)
//   http://ssl.ohmsha.co.jp/cgi-bin/menu.cgi?ISBN=4-274-07950-3

class PhotonFilter;

//

/**
 * This is the photon
 * The power is not compressed so the
 * size is 28bytes
 */
typedef struct Photon {
	float pos[3];					// photon position
	short plane;					// splitting plane for kd-tree
	unsigned char theta, phi;   	// incoming direction
	float power[3];					// photon power (uncompressed)
} Photon;

/**
 * This structure is used only to locate the
 * nearest photons
 */
typedef struct NearestPhotons {
	int max;
	int found;
	int got_heap;
	float pos[3];
    float normal[3];
	float *dist2;
	const Photon **index;
} NearestPhotons;


/**
 * This is the Photon_map class
 */
class Photon_map {
public:
	Photon_map(int max_phot);
	~Photon_map();

    void SetFilter(const PhotonFilter* pFilter);
    
	void store(
		const float power[3],		// photon power
		const float pos[3],  		// photon position
		const float dir[3]);		// photon direction

	void scale_photon_power(
		const float scale);			// 1/(number of emitted photons)

	void balance(void);				// balance the kd-tree (before use!)

	void irradiance_estimate(
		float irrad[3],				// returned irradiance
		const float pos[3],			// surface position
		const float normal[3],		// surface normal at pos
		const float max_dist,		// max distance to look for photons
		const int nphotons) const;	// number of photon to use
	
	void locate_photons(
		NearestPhotons* const np,	// np is used to locate the photons
        const int index) const;     // call with index = 1
	
	void photon_dir(
		float* dir,					// dir of photon (returned)
		const Photon* p) const;     // the photon
	
private:
	
	void balance_segment(
		Photon** pbal,
		Photon** porg,
		const int index,
		const int start,
		const int end);

	void median_split(
		Photon** p,
		const int start,
		const int end,
		const int median,
		const int axis);

	Photon* photons;

	int stored_photons;
	int half_stored_photons;
	int max_photons;
	int prev_scale;

	float costheta[256];
	float sintheta[256];
	float cosphi[256];
	float sinphi[256];

	float bbox_min[3];		// photons全体のbbox_min
	float bbox_max[3];		// photons全体のbbox_max
    
    const PhotonFilter* pFilter_;
};

