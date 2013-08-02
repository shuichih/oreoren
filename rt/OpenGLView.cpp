#ifdef __APPLE__

#include "OpenGLView.h"
#include "Config.h"
#include "BVH.h"
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#include <cstdio>

static OpenGLView* s_pView = 0;

//

static void Display()
{
    s_pView->DrawFrame();
}

static void Idle()
{
    glutPostRedisplay();
}

//

OpenGLView::OpenGLView(Config* pConfig, BVH* pBVH)
: initialized_(false)
, pConf_(pConfig)
, pBVH_(pBVH)
{
}

OpenGLView::~OpenGLView()
{
}

bool OpenGLView::Init(int w, int h)
{
    if (initialized_) {
        return true;
    }

    s_pView = this;

    width_ = w;
    height_ = h;
    int argc = 0;
    char** argv = NULL;
    glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(w, h);
	glutCreateWindow("Renderer");
	glutDisplayFunc(Display);
    //glutIdleFunc(Idle);

    initialized_ = true;
    
    return true;
}

bool OpenGLView::Present(u8* pColorBuf)
{
    if (!initialized_) {
        printf("The view is not initialized.\n");
        return false;
    }

    pColorBuf_ = pColorBuf;
    
	glutMainLoop();

    return true;
}

bool OpenGLView::DrawFrame()
{
    DrawToBuffer(pColorBuf_);
    DrawDebugStuff();

    glutSwapBuffers();
    
    return true;
}

void OpenGLView::DrawToBuffer(u8* pColorBuf)
{
    // Display
    int w = width_;
    int h = height_;
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

void OpenGLView::DrawDebugStuff()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = width_ / (float)height_;
    
    float fov = (float)atan2(0.5135, (double)pConf_->camera.direction.length());
    float fov_deg = Rad2Deg(fov);
    //fov_deg = 29;
    gluPerspective(fov_deg, aspect, 140.0, 1000.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    Vec3 pos = pConf_->camera.position;
    Vec3 at = pos + pConf_->camera.direction;
    
    // Vec3 up
    // const Vec3 proj_plane_axis_x = Vec3(fovY, 0.f, 0.f);
    // const Vec3 proj_plane_axis_y = (proj_plane_axis_x % camRay.d).normalize() * fovY;
    
    gluLookAt(pos.x, pos.y, pos.z, at.x, at.y, at.z, 0.0, 1.0, 0.0);

    glDisable(GL_LIGHTING);
    
    if (pConf_->bvhConf.draw) {
        DrawBVH(pBVH_, 0);
    }
    if (pConf_->drawBBox) {
        DrawBBox();
    }
    
    // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
}

void OpenGLView::DrawBVH(const IShape* pShape, int depth)
{
    if (pShape == NULL || depth >= pConf_->bvhConf.drawDepth) {
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

void OpenGLView::DrawBBox()
{
    for (int i=0; i<pConf_->scene.GetShapeNum(); i++) {
        const IShape* pShape = pConf_->scene.GetShape(i);
        if (!pShape->IsBVH()) {
            const BBox& bbox = pShape->BoundingBox();
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

#endif
