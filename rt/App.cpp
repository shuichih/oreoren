//
//  App.cpp
//  rt
//
//  Created by 秀一 林 on 11/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <GLUT/glut.h>
#include "App.h"
#include "smallpt_fmt.h"
#include <algorithm>
#include "PhotonMap.h"
#include "PhotonMapRenderer.h"

using namespace std;

static App* s_pApp = 0;

void DrawScene();
void Idle();

App::App()
{

}

void App::Run(int argc, const char * argv[], int w, int h, Mode mode)
{
    s_pApp = this;
    Init(argc, argv, w, h, mode);
}

void App::Init(int argc, const char * argv[], int w, int h, Mode mode)
{
    mode_ = mode;
    w_ = w;
    h_ = h;
    //pColorBuf_ = float[w_ * h_ * 3];
    //pPhotonMap_ = new Photon_map(500000);

    glutInit(&argc, (char**)argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(w, h);
	glutCreateWindow("glut window");
	glutDisplayFunc(DrawScene);
    //glutIdleFunc(Idle);
	glutMainLoop();
}

void App::Update()
{
    static unsigned char* pColorBuf = 0;

    //pColorBuf = rt(w_, h_, 2);
    Photon_map photonMap(500000);
    PhotonMapRenderer renderer;
    pColorBuf = renderer.Run(&photonMap, 500000, w_, h_);

#if 1
    glMatrixMode(GL_PROJECTION);
    glOrtho(0, w_, h_, 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
#else
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, 1.0, 0.1, 1000.0);
    glMatrixMode(GL_MODELVIEW);
#endif
    glDisable(GL_DITHER);

    glEnable(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
                
    glTexImage2D(GL_TEXTURE_2D, 0, 3, w_, h_, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, pColorBuf);
    
    //gentextureとかいらないのは何なのかね
    
    glClearColor(0.5, 0.5, 0.5, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glBegin(GL_POLYGON);
    glDisable(GL_CULL_FACE);
    glTexCoord2f(1.0, 1.0); // texcoordとvertexはこの順じゃないと駄目
    glVertex3f(w_, h_, 0.0);
    glTexCoord2f(0.0, 1.0);
    glVertex3f(0, h_, 0.0);
    glTexCoord2f(0.0, 0.0);
    glVertex3f(0, 0, 0.0);
    glTexCoord2f(1.0, 0.0);
    glVertex3f(w_, 0, 0.0);
    glEnd();
    
    glDisable(GL_TEXTURE_2D);

    glutSwapBuffers();
}

void DrawScene()
{
    if (s_pApp) {
        s_pApp->Update();
    }
}

void Idle()
{
    glutPostRedisplay();
}
