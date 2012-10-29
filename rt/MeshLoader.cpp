#include "MeshLoader.h"
#include "objLoader/objLoader.h"
#include <cstdio>
#include "Common.h"
#include "Scene.h"
#include <cassert>


ObjLoader::ObjLoader()
{
    pLoader_ = new objLoader();
}

ObjLoader::~ObjLoader()
{
    delete pLoader_;
}

void ObjLoader::SetFaceReverse(bool reverse)
{
    faceReverse_ = reverse;
}

Mesh* ObjLoader::Load(const char* pFilePath)
{
    // objをロードする
    if (!pLoader_->load(pFilePath)) {
        return NULL;
    }
    
    // 頂点数取得
    int nVertices = pLoader_->vertexCount;
    
    // 三角形換算の面数を算出
    int nFaces = 0;
    for (int i = 0; i < pLoader_->faceCount; i++) {
        obj_face* pObjFace = pLoader_->faceList[i];
        if (pObjFace->vertex_count > 5) {
            printf("Face[%d] has %d(>5) vertices\n", i, pObjFace->vertex_count);
            assert(false);
        }
        nFaces += (pObjFace->vertex_count - 2);
    }
    
    // 面と頂点の領域確保
    pMesh_ = new Mesh(nVertices, nFaces);
    
    // 頂点のコピー
    for (int i = 0; i < pLoader_->vertexCount; i++) {
        obj_vector* pVert = pLoader_->vertexList[i];
        pMesh_->pVertices[i].pos = Vec3((real)pVert->e[0], (real)pVert->e[1], (real)pVert->e[2]);
    }
    
    // 面法線と頂点法線を自前算出
    pVertNs_ = new Vec3[pLoader_->vertexCount];
    pFaceNs_ = new Vec3[pLoader_->faceCount];
    for (int i = 0; i < pLoader_->faceCount; i++)
    {
        obj_face* pObjFace = pLoader_->faceList[i];
        obj_vector* pVert0 = pLoader_->vertexList[pObjFace->vertex_index[0]];
        obj_vector* pVert1 = pLoader_->vertexList[pObjFace->vertex_index[1]];
        obj_vector* pVert2 = pLoader_->vertexList[pObjFace->vertex_index[2]];
        Vec3 p0((real)pVert0->e[0], (real)pVert0->e[1], (real)pVert0->e[2]);
        Vec3 p1((real)pVert1->e[0], (real)pVert1->e[1], (real)pVert1->e[2]);
        Vec3 p2((real)pVert2->e[0], (real)pVert2->e[1], (real)pVert2->e[2]);
        Vec3 s0 = p1 - p0;
        Vec3 s1 = p2 - p0;
        pFaceNs_[i] = (s0 % s1).normalize();
        if (faceReverse_) {
            pFaceNs_[i] *= -1;
        }
        
        for (int iVert = 0; iVert < pObjFace->vertex_count; iVert++) {
            int vidx = pObjFace->vertex_index[iVert];
            pVertNs_[vidx] += pFaceNs_[i];
        }
    }
    for (int i = 0; i < pLoader_->vertexCount; i++) {
        pVertNs_[i].normalize();
    }
    
    // 面を三角形に分割して格納
    u32 iFace = 0;
    for (int iObjFace=0; iObjFace<pLoader_->faceCount; iObjFace++) {
        obj_face* pObjFace = pLoader_->faceList[iObjFace];
        Vec3& rFaceNorm = pFaceNs_[iObjFace];
        
        for (int iObjVert=0; iObjVert<pObjFace->vertex_count-2; iObjVert++) {
            ProcessFace(iFace++, pObjFace, rFaceNorm, 0, 1+iObjVert, 2+iObjVert);
        }
    }
    
    printf("%d faces, %d vertices.\n", pMesh_->nFaces, pMesh_->nVertices);
    
    return pMesh_;
}

void ObjLoader::ProcessFace(
    int iFace,
    obj_face* pObjFace,
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
        int vidx = pObjFace->vertex_index[vertLocalIdx];
        pMesh_->pFaces[iFace].indices[iVert] = vidx;
        Vertex* pV = &pMesh_->pVertices[vidx];
        
        int nidx = pObjFace->normal_index[vertLocalIdx];
        if (nidx != -1 && nidx < pLoader_->normalCount) {
            obj_vector* pNorm = pLoader_->normalList[nidx];
            pV->normal.x = (float)pNorm->e[0];
            pV->normal.y = (float)pNorm->e[1];
            pV->normal.z = (float)pNorm->e[2];
            pV->normal.normalize();
            if (faceReverse_) {
                pV->normal *= -1;
            }
        } else {
            pV->normal = pVertNs_[vidx];
        }
    }
}
