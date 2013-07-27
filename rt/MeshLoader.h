#ifndef MeshLoader_h
#define MeshLoader_h

#include "Common.h"

class Mesh;
class objLoader;
struct obj_face;

class ObjLoader
{
public:
    ObjLoader();
    ~ObjLoader();
    
    Mesh* Load(const char* pFilePath);
    void ProcessFace(int iFace, obj_face* pObjFace, Vec3& rObjFaceNorm, int i0, int i1, int i2);
    void SetFaceReverse(bool reverse);
    
private:
    bool faceReverse_;
    objLoader* pLoader_;
    
    // work
    Mesh* pMesh_;
    Vec3* pVertNs_;
    Vec3* pFaceNs_;
};

#endif
