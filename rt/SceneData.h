#ifndef _SceneData_H_
#define _SceneData_H_

#include "Scene.h"

extern Sphere g_spheres[];
extern int g_nSpheres;


// 座標系は右手系, X軸:右方向 Y軸:上方向 Z軸:手前方向
Sphere g_spheres[] = {
    //Scene: radius, position, emission, color, material
    Sphere(1e5, Vec( 1e5+1,  40.8,      81.6),       Vec(),Vec(.75,.25,.25),DIFF),    //Left
    Sphere(1e5, Vec(-1e5+99, 40.8,      81.6),       Vec(),Vec(.25,.25,.75),DIFF),    //Rght
    Sphere(1e5, Vec(50,      40.8,      1e5),        Vec(),Vec(.75,.75,.75),DIFF),    //Back
    Sphere(1e5, Vec(50,      40.8,      -1e5+170),   Vec(),Vec(),           DIFF),    //Frnt
    Sphere(1e5, Vec(50,      1e5,       81.6),       Vec(),Vec(.75,.75,.75),DIFF),    //Botm
    Sphere(1e5, Vec(50,      -1e5+81.6, 81.6),       Vec(), Vec(.75,.75,.75),DIFF),   //Top
    Sphere(16.5,Vec(27,      16.5,      47),         Vec(),Vec(1,1,1)*.999, SPEC),    //Mirr
    Sphere(16.5,Vec(73,      16.5,      78),         Vec(),Vec(1,1,1)*.999, REFR),    //Glas
  //Sphere(600, Vec(50,      681.6-.27, 81.6),       Vec(12,12,12), Vec(),  DIFF)     //Lite
};


int g_nSpheres = sizeof(g_spheres) / sizeof(g_spheres[0]);

#endif
