//
//  App.cpp
//  rt
//
//  Created by 秀一 林 on 11/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#include <algorithm>
#include <cassert>
#include "App.h"
#include "smallpt_fmt.h"
#include "PhotonMap.h"
#include "PhotonMapRenderer.h"
#include "PhotonMapRenderer2.h"
#include "PostEffect.h"
#include "Scene.h"
#include "BVH.h"
#include "QBVH.h"
#include "RayTracingRenderer.h"


using namespace std;

static App* s_pApp = 0;

void Display();
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
: pBVH_(NULL)
, pQBVH_(NULL)
{

}

App::~App()
{
    delete pBVH_;
    delete pQBVH_;
}


void App::Run(int argc, const char * argv[])
{
    s_pApp = this;
    
    Init(argc, argv);
    
	glutMainLoop();
}

void App::Init(int argc, const char * argv[])
{
    // Configロード
    //const char* path = (argc >= 2) ? argv[1] : "~/Dev/rt/rt/SceneFiles/cornell_box.scene";
    //const char* path = (argc >= 2) ? argv[1] : "~/Dev/rt/rt/SceneFiles/venus.scene";
    const char* path = (argc >= 2) ? argv[1] : "~/Dev/rt/rt/SceneFiles/simple.scene";
    config.Load(path);
    printf("--------------------------------\n\n");
    
    // GL初期化
    glutInit(&argc, (char**)argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(config.windowWidth, config.windowHeight);
	glutCreateWindow("Renderer");
	glutDisplayFunc(Display);
    //glutIdleFunc(Idle);
    
    // レンダラにConfig内容を設定
    if (config.rendererType == RTYPE_SIMPLE_RT) {
        pRenderer_ = new RayTracingRenderer();
    }
    else if (config.rendererType == RTYPE_PHOTON_MAP) {
        pRenderer_ = new PhotonMapRenderer();
    }
    else if (config.rendererType == RTYPE_PHOTON_MAP2) {
        pRenderer_ = new PhotonMapRenderer2();
    }
    pRenderer_->SetConfig(config);
}

void App::BuildBVH()
{
    if (config.bvhConf.build) {
        Timer timer("Building BVH");
        config.scene.BuildBVH(config.bvhConf.type);
    }
}

void App::ConvertToUint(u8* pColorBuf, Vec3* pRealColorBuf)
{
    // Convert to u8 format
    int w = config.windowWidth;
    int h = config.windowHeight;
    for (int i = 0, j=0; i < (w*h); ++i, j+=4)
    {
        pColorBuf[j+0] = (u8)(toInt(pRealColorBuf[i].x)); // including Gamma Correction
        pColorBuf[j+1] = (u8)(toInt(pRealColorBuf[i].y));
        pColorBuf[j+2] = (u8)(toInt(pRealColorBuf[i].z));
        pColorBuf[j+3] = 255;
    }
    
}

void App::Render()
{
    static bool first = true;
    
    if (!first) return;
    
    int w = config.windowWidth;
    int h = config.windowHeight;
    
    Vec3* pRealColorBuf = new Vec3[w * h];
    u8* pColorBuf = new u8[w * h * sizeof(char)*4];
    
    {
        Timer timer;
        
        BuildBVH();
       
        RenderScene(pRealColorBuf);
        
        timer.PrintElapsed("Total rendering time: ");
    }
    
    ConvertToUint(pColorBuf, pRealColorBuf);
    DrawToBuffer(pColorBuf);
    DrawDebugStuff();
    
    glutSwapBuffers();
    
    delete [] pRealColorBuf;
    delete [] pColorBuf;
    
    first = false;
}

void App::RenderScene(Vec3* pRealColorBuf)
{
    // Render using photonmap
    pRenderer_->Run(pRealColorBuf, config.scene);
    
    // Tone Mapping
    if (config.postEffect.toneMapEnabled) {
        ToneMap toneMap;
        toneMap.SetKeyValue(config.postEffect.toneMapKeyValue);
        toneMap.Apply(pRealColorBuf, config.windowWidth, config.windowHeight);
    }
}

void App::DrawToBuffer(u8* pColorBuf)
{
    int w = config.windowWidth;
    int h = config.windowHeight;
    
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

}

void App::DrawDebugStuff()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = config.windowWidth / (float)config.windowHeight;
    
    float fov = (float)atan2(0.5135, (double)config.camera.direction.length());
    float fov_deg = Rad2Deg(fov);
    //fov_deg = 29;
    gluPerspective(fov_deg, aspect, 140.0, 1000.0);
//    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Vec3 pos = config.camera.position;
    Vec3 at = pos + config.camera.direction;
    
   // Vec3 up
   // const Vec3 proj_plane_axis_x = Vec3(fovY, 0.f, 0.f);
   // const Vec3 proj_plane_axis_y = (proj_plane_axis_x % camRay.d).normalize() * fovY;
    
    gluLookAt(pos.x, pos.y, pos.z, at.x, at.y, at.z, 0.0, 1.0, 0.0);

    glDisable(GL_LIGHTING);
    
    if (config.bvhConf.draw) {
        DrawBVH(pBVH_, 0);
    }
    if (config.drawBBox) {
        DrawBBox();
    }
    
   // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    
}

void App::DrawBVH(const IShape* pShape, int depth)
{
    if (pShape == NULL || depth >= config.bvhConf.drawDepth) {
        return;
    }
    
    GLfloat color[4];
    if (pShape->IsBVH()) {
        color[0] = 0.f; color[1] = 0.f; color[2] = 1.f; color[3] = 1.f;
    } else {
        color[0] = 0.f; color[1] = 1.f; color[2] = 0.f; color[3] = 1.f;
    }
    const BBox& bbox = pShape->BoundingBox();
    Vec3 center = bbox.Center();
    Vec3 size = bbox.Size();
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glColor3f(color[0], color[1], color[2]);
    glTranslatef(center.x, center.y, center.z);
    glScalef(size.x, size.y, size.z);
    glutWireCube(1);
    glPopMatrix();
    
    if (pShape->IsBVH()) {
        BVH* pBVH = (BVH*)pShape;
        DrawBVH(pBVH->pLeft_, depth+1);
        DrawBVH(pBVH->pRight_, depth+1);
    }
}

void App::DrawBBox()
{
    std::vector<const IShape*>& shapes = config.scene.shapes_;
    for (int i=0; i<shapes.size(); i++) {
        if (!shapes[i]->IsBVH()) {
            const BBox& bbox = shapes[i]->BoundingBox();
            Vec3 center = bbox.Center();
            Vec3 size = bbox.Size();
            
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glColor3f(1, 0, 0);
            glTranslatef(center.x, center.y, center.z);
            glScalef(size.x, size.y, size.z);
            glutWireCube(1);
            glPopMatrix();
        }
    }
}

void Display()
{
    s_pApp->Render();
}

void Idle()
{
    glutPostRedisplay();
}
