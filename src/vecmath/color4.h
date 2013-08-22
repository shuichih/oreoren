#ifndef SEVERE3D_COLOR4_H
#define SEVERE3D_COLOR4_H

//#include "severe3d/base.h"
#include "vecmath_def.h"

/**
 * ARGB色を表現するクラスです。
 */
class Color4
{
public:

	/** 透明度 (透明 < 不透明) */
	float a;

	/** 赤成分 */
	float r;

	/** 緑成分 */
	float g;

	/** 青成分 */
	float b;

	Color4(): a(0), r(0), g(0), b(0) {}
	Color4(float sr, float sg, float sb): a(1.0f), r(sr), g(sg), b(sb) {}
	Color4(float sa, float sr, float sg, float sb): a(sa), r(sr), g(sg), b(sb) {}

	inline uint32 to32BitARGB() const {
		// TODO
		uint32 uc_a = static_cast<uint32>(a * 255);
		uint32 uc_r = static_cast<uint32>(r * 255);
		uint32 uc_g = static_cast<uint32>(g * 255);
		uint32 uc_b = static_cast<uint32>(b * 255);

		return (uc_a << 24) | (uc_r << 16) | (uc_g << 8) | uc_b;
	}

	inline void from32BitARGB(int32 color32) {
		a = ((color32 >> 24) & 0x000000ff) / 255.0f;
		r = ((color32 >> 16) & 0x000000ff) / 255.0f;
		g = ((color32 >> 8)  & 0x000000ff) / 255.0f;
		b = ( color32        & 0x000000ff) / 255.0f;
	}

	inline void mul(const Color4& color) {
		r *= color.r;
		g *= color.g;
		b *= color.b;
	}

	inline void mul(const Color4& c1, const Color4& c2) {
		r = c1.r * c2.r;
		g = c1.g * c2.g;
		b = c1.b * c2.b;
	}

    inline void scaleAdd(const Color4& c1, const Color4& c2, float scale) {
		// alias safe
        r = c1.r + scale * c2.r;
        g = c1.g + scale * c2.g;
        b = c1.b + scale * c2.b;
    }

    inline void scaleAdd(const Color4& c1, float scale) {
        r += scale * c1.r;
        g += scale * c1.g;
        b += scale * c1.b;
    }

	inline void scale(float scale) {
        r *= scale;
        g *= scale;
        b *= scale;
	}

	inline void add(const Color4& c1) {
		r += c1.r;
		g += c1.g;
		b += c1.b;
	}

	inline void sub(const Color4& c1) {
		r -= c1.r;
		g -= c1.g;
		b -= c1.b;
	}

	inline Color4& clamp() {
		if (r < 0) r = 0; else if (r > 1.0f) r = 1.0f;
		if (g < 0) g = 0; else if (g > 1.0f) g = 1.0f;
		if (b < 0) b = 0; else if (b > 1.0f) b = 1.0f;
		return *this;
	}

	inline void zero() {
		a = 0;
		r = 0;
		g = 0;
		b = 0;
	}

	// NOTE:コピーコンストラクタおよび = 演算子はコンパイラによって生成されます。

	inline unsigned char a_byte() {
		return static_cast<unsigned char>(a * 255);
	}
	inline unsigned char r_byte() {
		return static_cast<unsigned char>(r * 255);
	}
	inline unsigned char g_byte() {
		return static_cast<unsigned char>(g * 255);
	}
	inline unsigned char b_byte() {
		return static_cast<unsigned char>(b * 255);
	}

    Color4& operator+=(const Color4& c1) {
        add(c1);
        return *this;
    }
    Color4& operator-=(const Color4& c1) {
        sub(c1);
        return *this;
    }
    Color4& operator*=(const Color4& c1) {
        mul(c1);
        return *this;
    }
    Color4& operator*=(float s) {
        scale(s);
        return *this;
    }
    Color4 operator+(const Color4& c1) const {
        return (Color4(*this)).operator+=(c1);
    }
    Color4 operator-(const Color4& c1) const {
        return (Color4(*this)).operator-=(c1);
    }
    Color4 operator*(const Color4& c1) const {
        return (Color4(*this)).operator*=(c1);
    }
    Color4 operator*(float s) const {
        return (Color4(*this)).operator*=(s);
    }

};


// グローバル演算子
template <class T>
inline
severe3d::Color4 operator*(T s, const severe3d::Color4& c1) {
    return c1 * s;
}


#endif // SEVERE3D_COLOR4_H
