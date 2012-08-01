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
    Sphere(1e5f, Vec( 1+1e5f,  40.8f,      81.6f),     Vec(.75f, .25f, .25f),    DIFF),    //Left
    Sphere(1e5f, Vec(99-1e5f,  40.8f,      81.6f),     Vec(.25f, .25f, .75f),    DIFF),    //Rght
    Sphere(1e5f, Vec(50.f,     40.8f,      1e5f),      Vec(.75f, .75f, .75f),    DIFF),    //Back
    //Sphere(1e5, Vec(50f,      40.8f,      170-1e5f), Vec(),                     DIFF),    //Frnt
    Sphere(1e5f, Vec(50.f,     40.8f,      170-1e5f),  Vec(.75f, .75f, .75f),    DIFF),    //Frnt
    Sphere(1e5f, Vec(50.f,      1e5f,      81.6f),     Vec(.75f, .75f, .75f),    DIFF),    //Botm
    Sphere(1e5f, Vec(50.f,     81.6f-1e5f, 81.6f),     Vec(.75f, .75f, .75f),    DIFF),    //Top
    Sphere(16.5f, Vec(27.f,    16.5f,      47.f),      Vec(1.f, 1.f, 1.f)*.999f, SPEC),    //Mirr
    Sphere(16.5f,Vec(73.f,     16.5f,      78.f),      Vec(1.f, 1.f, 1.f)*.999f, REFR),    //Glas
};

Triangle g_triangles[] = {
//    Triangle(Vec(50, 40.8+20, 85), Vec(50-20, 40.8-20, 85), Vec(50+20, 40.8-20, 85), Vec(.80, .80, .20), DIFF),
//Triangle(Vec(50, 40.8+20, 85), Vec(50+20, 40.8-20, 85), Vec(50-20, 40.8-20, 85), Vec(.80, .80, .20), DIFF),
};

int g_nSpheres = sizeof(g_spheres) / sizeof(g_spheres[0]);
int g_nTriangles = sizeof(g_triangles) / sizeof(g_triangles[0]);

Shape** g_shapes = NULL;
int g_nShapes = 0;

#endif
