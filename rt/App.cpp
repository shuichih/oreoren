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
#include "PostEffect.h"
#include "MeshLoader.h"
#include "Scene.h"


using namespace std;

static App* s_pApp = 0;

void DrawScene();
void Idle();

//

inline real clamp(real x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(real x)
{
#ifdef USE_FLOAT
    return int(powf(clamp(x), 1/2.2) * 255 + .5);
#else
    return int(pow(clamp(x), 1/2.2) * 255 + .5);
#endif
}

//

App::App()
{

}

App::~App()
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

    glutInit(&argc, (char**)argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(w, h);
	glutCreateWindow("glut window");
	glutDisplayFunc(DrawScene);
    //glutIdleFunc(Idle);
    
    // objロード
    /*
    const char* OBJ_FILE = "/Users/shuichih/Dev/rt/venusm.obj";
    ObjLoader loader;
    Mesh* pMesh = loader.Load(OBJ_FILE);
	if (pMesh == NULL)
	{
        printf("File Load Error: %s\n", OBJ_FILE);
		return;
	}
    pMesh->scale(0.01, 0.01, 0.01);
    */
    
	glutMainLoop();
}

void App::Update()
{
    static u8* pColorBuf = 0;

    // Render using photonmap
    PhotonMapRenderer renderer;
    PhotonMapRenderer::Config config = renderer.GetDefaultConfig();
    config.screenWidth = w_;
    config.screenHeight = h_;
    config.nPhotons = 100000;
    config.nEstimatePhotons = 200;
    config.estimateDist = 10.f;
    config.nSubPixelsSqrt = 1;
    renderer.SetConfig(config);
    
    time_t startTime, endTime;
    time(&startTime);
    
    Vec* pRealColorBuf = renderer.Run();
    
    // Tone Mapping
    ToneMap toneMap;
    toneMap.SetKeyValue(0.045);
    //toneMap.SetDelta(0.01);
    //toneMap.Apply(pRealColorBuf, w_, h_);
    
    // Convert to u8 format
    delete pColorBuf;
    pColorBuf = new u8[w_ * h_ * sizeof(char)*4];
    for (int i = 0, j=0; i < (w_*h_); ++i, j+=4)
    {
        pColorBuf[j+0] = (u8)(toInt(pRealColorBuf[i].x)); // including Gamma Correction
        pColorBuf[j+1] = (u8)(toInt(pRealColorBuf[i].y));
        pColorBuf[j+2] = (u8)(toInt(pRealColorBuf[i].z));
        pColorBuf[j+3] = 255;
    }
    delete [] pRealColorBuf;
    
    time(&endTime);
    u32 elapsed = (u32)difftime(endTime, startTime);
    printf("rendering time = %dm %ds\n", elapsed/60, elapsed%60);

    // Display
    glMatrixMode(GL_PROJECTION);
    glOrtho(0, w_, h_, 0, -1, 1);
    //gluPerspective(60.0, 1.0, 0.1, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glDisable(GL_DITHER);

    glEnable(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
                
    glTexImage2D(GL_TEXTURE_2D, 0, 3, w_, h_, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, pColorBuf);
    
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
