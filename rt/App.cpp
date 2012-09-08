//
//  App.cpp
//  rt
//
//  Created by 秀一 林 on 11/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <GLUT/glut.h>
#include <algorithm>
#include "App.h"
#include "smallpt_fmt.h"
#include "PhotonMap.h"
#include "PhotonMapRenderer.h"
#include "PostEffect.h"
#include "Scene.h"


using namespace std;

static App* s_pApp = 0;

void DrawScene();
void Idle();

//

inline real clamp(real x)
{
    return x < 0.f ? 0.f : x > 1.f ? 1.f : x;
}

inline int toInt(real x)
{
#ifdef USE_FLOAT
    return int(powf(clamp(x), 1/2.2f) * 255.f + .5f);
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


void App::Run(int argc, const char * argv[])
{
    s_pApp = this;
    
    Init(argc, argv);
    
	glutMainLoop();
}

void App::Init(int argc, const char * argv[])
{
    const char* path = (argc >= 2) ? argv[1] : "~/Dev/rt/rt/SceneFiles/cornell_box.scene";
    config.Load(path);
    
    glutInit(&argc, (char**)argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(config.windowWidth, config.windowHeight);
	glutCreateWindow("Renderer");
	glutDisplayFunc(DrawScene);
    //glutIdleFunc(Idle);
    
    PhotonMapRenderer::Config pmconfig = renderer_.GetDefaultConfig();
    config.ApplyTo(pmconfig);
    renderer_.SetConfig(pmconfig);
    Ray cam(config.camera.position, config.camera.direction);
    renderer_.SetCamera(cam, config.camera.fovY);
}
    

void App::Update()
{
    static u8* pColorBuf = 0;

    int w = config.windowWidth;
    int h = config.windowHeight;
    
    // Render using photonmap
    
    time_t startTime, endTime;
    time(&startTime);
    
    Vec3* pRealColorBuf = renderer_.Run(config.scene);
    
    // Tone Mapping
    if (config.postEffect.toneMapEnabled) {
        ToneMap toneMap;
        toneMap.SetKeyValue(config.postEffect.toneMapKeyValue);
        toneMap.Apply(pRealColorBuf, w, h);
    }
    
    // Convert to u8 format
    delete pColorBuf;
    pColorBuf = new u8[w * h * sizeof(char)*4];
    for (int i = 0, j=0; i < (w*h); ++i, j+=4)
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
    glOrtho(0, w, h, 0, -1, 1);
    //gluPerspective(60.0f, 1.0f, 0.1f, 1000.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glDisable(GL_DITHER);

    glEnable(GL_TEXTURE_2D);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
                
    glTexImage2D(GL_TEXTURE_2D, 0, 3, w, h, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, pColorBuf);
    
    glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glBegin(GL_POLYGON);
    glDisable(GL_CULL_FACE);
    glTexCoord2f(1.0f, 1.0f); // texcoordとvertexはこの順じゃないと駄目
    glVertex3f(w, h, 0.0f);
    glTexCoord2f(0.0f, 1.0f);
    glVertex3f(0.f, h, 0.0f);
    glTexCoord2f(0.0f, 0.0f);
    glVertex3f(0.f, 0.f, 0.f);
    glTexCoord2f(1.0f, 0.0f);
    glVertex3f(w, 0.f, 0.0f);
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
