#ifndef _Common_H_
#define _Common_H_

struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
    
    union
    {
        struct
        {
            double x, y, z;
        };
        double e[3];
    };
    // position, also color (r,g,b)
    Vec(double x_=0, double y_=0, double z_=0){
        x=x_;
        y=y_;
        z=z_;
    }
    Vec operator+(const Vec &b) const {
        return Vec(x+b.x,y+b.y,z+b.z);
    }
    Vec operator-(const Vec &b) const {
        return Vec(x-b.x,y-b.y,z-b.z);
    }
    Vec operator*(double b) const {
        return Vec(x*b,y*b,z*b);
    }
    Vec mult(const Vec &b) const {
        return Vec(x*b.x,y*b.y,z*b.z);
    }
    Vec& norm(){
        return *this = *this * (1/sqrt(x*x+y*y+z*z));
    }
    double dot(const Vec &b) const {
        return x*b.x+y*b.y+z*b.z;
    }
    // cross:
    Vec operator%(Vec&b){
        return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
    }
};

struct Ray {
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};


#endif

