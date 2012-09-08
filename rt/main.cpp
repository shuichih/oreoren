//
//  main.cpp
//  rt
//
//  Created by 秀一 林 on 11/28/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
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

