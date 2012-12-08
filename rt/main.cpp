//
//  main.cpp
//  rt
//
//  Copyright 2012 Suichi Hayashi. All rights reserved.
//
#include <iostream>
#include "App.h"

#include <GLUT/glut.h>

void drawScene();

int main(int argc, const char * argv[])
{
    //printf("%lu %lu %lu %lu", sizeof(int), sizeof(long), sizeof(long long), sizeof(void*));

    App* pApp = new App();
    
    pApp->Run(argc, argv);

    delete pApp;
    
	return 0;
}


