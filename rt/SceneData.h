#ifndef SceneData_H
#define SceneData_H

#include "Scene.h"


extern Shape** g_shapes;
extern int g_nShapes;

extern Sphere g_spheres[];
extern int g_nSpheres;

// 座標系は右手系, X軸:右方向 Y軸:上方向 Z軸:手前方向
Sphere g_spheres[] = {
    //Scene: radius, position, emission, color, material
    Sphere(1e5, Vec( 1+1e5,  40.8,      81.6),    Vec(.75,.25,.25),DIFF),    //Left
    Sphere(1e5, Vec(99-1e5,  40.8,      81.6),    Vec(.25,.25,.75),DIFF),    //Rght
    Sphere(1e5, Vec(50,      40.8,      1e5),     Vec(.75,.75,.75),DIFF),    //Back
    //Sphere(1e5, Vec(50,      40.8,      170-1e5), Vec(),           DIFF),    //Frnt
    Sphere(1e5, Vec(50,      40.8,      170-1e5), Vec(.75,.75,.75),           DIFF),    //Frnt
    Sphere(1e5, Vec(50,      1e5,       81.6),    Vec(.75,.75,.75),DIFF),    //Botm
    Sphere(1e5, Vec(50,      81.6-1e5,  81.6),    Vec(.75,.75,.75),DIFF),    //Top
//    Sphere(16.5,Vec(27,      16.5,      47),      Vec(1,1,1)*.999, SPEC),    //Mirr
//    Sphere(16.5,Vec(73,      16.5,      78),      Vec(1,1,1)*.999, REFR),    //Glas
};

Triangle g_triangles[] = {
//    Triangle(Vec(50, 40.8+20, 85), Vec(50-20, 40.8-20, 85), Vec(50+20, 40.8-20, 85), Vec(.80, .80, .20), DIFF),
    Triangle(Vec(50, 40.8+20, 85), Vec(50+20, 40.8-20, 85), Vec(50-20, 40.8-20, 85), Vec(.80, .80, .20), DIFF),
};

int g_nSpheres = sizeof(g_spheres) / sizeof(g_spheres[0]);
int g_nTriangles = sizeof(g_triangles) / sizeof(g_triangles[0]);

Shape** g_shapes = NULL;
int g_nShapes = 0;

#endif
