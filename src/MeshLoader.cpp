#include "MeshLoader.h"
#include "Common.h"
#include "ObjLoader.h"
#include "Scene.h"
#include <cstdio>
#include <cassert>

using namespace std;

MeshLoader::MeshLoader()
{
    pLoader_ = new ObjLoader();
}

MeshLoader::~MeshLoader()
{
    delete pLoader_;
}

void MeshLoader::SetFaceReverse(bool reverse)
{
    faceReverse_ = reverse;
}

Mesh* MeshLoader::Load(const char* pFilePath)
{
    // objをロードする
    if (!pLoader_->Load(pFilePath)) {
        return NULL;
    }
    
    // 頂点数取得
    int nObjVertices = pLoader_->NumVertices();
    int nObjFaces = pLoader_->NumFaces();
    const vector<Vec3>& rObjVertices = pLoader_->Vertices();
    const vector<ObjFace>& rObjFaces = pLoader_->Faces();
    
    // 三角形換算の面数を算出
    int nFaces = 0;
    for (int i = 0; i < nObjFaces; i++) {
        const ObjFace& rObjFace = rObjFaces[i];
        if (rObjFace.nVertices > 5) {
            printf("Face[%d] has %d(>5) vertices\n", i, rObjFace.nVertices);
            assert(false);
        }
        nFaces += (rObjFace.nVertices - 2);
    }
    
    // 面と頂点の領域確保
    pMesh_ = new Mesh(nObjVertices, nFaces, NULL);
    
    // 頂点のコピー
    for (int i = 0; i < nObjVertices; i++) {
        pMesh_->pVertices[i].pos = rObjVertices[i];
    }
    
    // 面法線と頂点法線を自前算出
    pVertNs_ = new Vec3[nObjVertices];
    pFaceNs_ = new Vec3[nObjFaces];
    for (int i = 0; i < nObjFaces; i++)
    {
        const ObjFace& rObjFace = rObjFaces[i];
        Vec3 p0 = rObjVertices[rObjFace.iVertex[0]];
        Vec3 p1 = rObjVertices[rObjFace.iVertex[1]];
        Vec3 p2 = rObjVertices[rObjFace.iVertex[2]];
        Vec3 s0 = p1 - p0;
        Vec3 s1 = p2 - p0;
        pFaceNs_[i] = (s0 ^ s1).normalize();
        if (faceReverse_) {
            pFaceNs_[i] *= -1;
        }
        
        for (int iVert = 0; iVert < rObjFace.nVertices; iVert++) {
            int vidx = rObjFace.iVertex[iVert];
            pVertNs_[vidx] += pFaceNs_[i];
        }
    }
    for (int i = 0; i < nObjVertices; i++) {
        pVertNs_[i].normalize();
    }
    
    // 面を三角形に分割して格納
    u32 iFace = 0;
    for (int iObjFace=0; iObjFace<nObjFaces; iObjFace++) {
        const ObjFace& rObjFace = rObjFaces[iObjFace];
        Vec3& rFaceNorm = pFaceNs_[iObjFace];
        
        for (int iObjVert=0; iObjVert<rObjFace.nVertices-2; iObjVert++) {
            ProcessFace(iFace++, rObjFace, rFaceNorm, 0, 1+iObjVert, 2+iObjVert);
        }
    }
    
    printf("%d faces, %d vertices.\n", pMesh_->nFaces, pMesh_->nVertices);
    
    return pMesh_;
}

void MeshLoader::ProcessFace(
    int iFace,
    const ObjFace& rObjFace,
    Vec3& rObjFaceNorm,
    int i0,
    int i1,
    int i2)
{
    pMesh_->pFaces[iFace].normal = rObjFaceNorm;
    
    int surfLocalIndices[3] = { i0, i1, i2 };
    if (faceReverse_) {
        surfLocalIndices[1] = i2;
        surfLocalIndices[2] = i1;
    }
    
    for (int iVert = 0; iVert < 3; iVert++) {
        int vertLocalIdx = surfLocalIndices[iVert];
        int vidx = rObjFace.iVertex[vertLocalIdx];
        pMesh_->pFaces[iFace].indices[iVert] = vidx;
        Vertex* pV = &pMesh_->pVertices[vidx];
        
        int nidx = rObjFace.iNormal[vertLocalIdx];
        if (nidx != -1 && nidx < pLoader_->NumNormals()) {
            pV->normal = pLoader_->Normals()[nidx];
            pV->normal.normalize();
            if (faceReverse_) {
                pV->normal *= -1;
            }
        } else {
            pV->normal = pVertNs_[vidx];
        }
    }
}
