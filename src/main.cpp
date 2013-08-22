//
//  main.cpp
//  rt
//
//  Copyright 2012 Suichi Hayashi. All rights reserved.
//
#include <iostream>
#include "App.h"

//#include <GLUT/glut.h>
#include "vecmath/tuple3.h"
#include "vecmath/vector3.h"

void printVec(const char* pStr, Vector3f& a);
void TestVector();

int main(int argc, const char * argv[])
{
    //printf("%lu %lu %lu %lu", sizeof(int), sizeof(long), sizeof(long long), sizeof(void*));

#if 1
    App* pApp = new App();
    
    pApp->Run(argc, argv);

    delete pApp;
#else
    TestVector();
#endif
    
	return 0;
}

void printVec(const char* pStr, Vector3f& a)
{
    printf("%8s %.8f %.8f %.8f\n", pStr, a.x, a.y, a.z);
}

void TestVector()
{
    Vector3f v[21];
    for (int i=0; i<sizeof(v)/sizeof(v[0]); i++) {
        v[i] = Vector3f((float)i, i*2+0.1f, i*3+0.2f);
    }
    v[0] = Vector3f(1, 2, 3);
    v[1].set(4, 5, 6);
    printVec("init", v[0]);
    printVec("set", v[1]);
    
    printVec("sub0", v[2]);
    printVec("sub1", v[3]);
    v[2].sub(v[3]);
    printVec("subA", v[2]);
    
    printVec("add0", v[4]);
    printVec("add1", v[5]);
    v[4].add(v[5]);
    printVec("addA", v[4]);
    
    printVec("scale", v[6]);
    v[6].scale(1.5f);
    printVec("scaleA", v[6]);
    
    printVec("negate", v[7]);
    v[7].negate();
    printVec("nagateA", v[7]);
    
    printVec("interp0", v[8]);
    printVec("interp1", v[9]);
    v[8].interpolate(v[8], v[9], 0.1f);
    printVec("interpA", v[8]);
    
    printVec("assignD", v[10]);
    printVec("assignS", v[11]);
    v[10] = v[11];
    printVec("assignA", v[10]);
    
    printVec("lengthSq", v[12]);
    float lengthSq = v[12].lengthSquared();
    printVec("lengthSq", v[12]);
    printf("lengthSq %f\n", lengthSq);
    
    printVec("length", v[13]);
    float length = v[13].length();
    printVec("length", v[13]);
    printf("length %f\n", length);
    
    v[15] = v[14];
    printVec("normal", v[14]);
    v[14].normalize();
    printVec("normal", v[14]);
    
    printVec("normalF", v[15]);
    v[15].normalizeFast();
    printVec("normalF", v[15]);
    
    printVec("dot", v[16]);
    printVec("dot", v[17]);
    float dotA = v[16].dot(v[17]);
    printf("dotA %f\n", dotA);
    
    printVec("cross", v[18]);
    printVec("cross", v[19]);
    v[20].cross(v[18], v[19]);
    printVec("crossA", v[20]);
    
    Vector3f crossA = v[18].cross(v[19]);
    printVec("crossA", crossA);
}
