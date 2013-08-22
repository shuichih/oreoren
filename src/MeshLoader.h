#ifndef MeshLoader_h
#define MeshLoader_h

#include "Common.h"

class Mesh;
class ObjLoader;
struct ObjFace;

class MeshLoader
{
public:
    MeshLoader();
    ~MeshLoader();
    
    Mesh* Load(const char* pFilePath);
    void ProcessFace(int iFace, const ObjFace& rObjFace, Vec3& rObjFaceNorm, int i0, int i1, int i2);
    void SetFaceReverse(bool reverse);
    
private:
    bool faceReverse_;
    ObjLoader* pLoader_;
    
    // work
    Mesh* pMesh_;
    Vec3* pVertNs_;
    Vec3* pFaceNs_;
};

#endif
