#include "MeshLoader.h"
#include "objLoader/objLoader.h"
#include <cstdio>
#include "Common.h"
#include "Scene.h"
#include <cassert>


ObjLoader::ObjLoader()
{
}

ObjLoader::~ObjLoader()
{
}

Mesh* ObjLoader::Load(const char* pFilePath)
{
    objLoader* pLoader = new objLoader();
    if (!pLoader->load(pFilePath)) {
        return NULL;
    }
    
    int nVertices = pLoader->vertexCount;
    
    // 三角形換算で面数を数える
    int nFaces = 0;
    for (int i = 0; i < pLoader->faceCount; i++) {
        obj_face* pFace = pLoader->faceList[i];
        //DebugPrint( L"f=%d¥n", i );
        if (pFace->vertex_count == 3) {
            nFaces++;
        } else if (pFace->vertex_count == 4) {
            nFaces += 2;
        } else {
            printf("Face[%d] has >4 vertices\n", i);
            assert(false);
        }
    }
    
    Mesh* pMesh = new Mesh(nVertices, nFaces);
    
    for (int i = 0; i < pLoader->vertexCount; i++) {
        obj_vector* pVert = pLoader->vertexList[i];
        pMesh->pVertices[i].pos = Vec((real)pVert->e[0], (real)pVert->e[1], (real)pVert->e[2]);
    }
    
    Vec* pVertNs = new Vec[pLoader->vertexCount];
    Vec* pFaceNs = new Vec[pLoader->faceCount];
    for (int i = 0; i < pLoader->faceCount; i++)
    {
        obj_face* pFace = pLoader->faceList[i];
        obj_vector* pVert0 = pLoader->vertexList[pFace->vertex_index[0]];
        obj_vector* pVert1 = pLoader->vertexList[pFace->vertex_index[1]];
        obj_vector* pVert2 = pLoader->vertexList[pFace->vertex_index[2]];
        Vec p0((real)pVert0->e[0], (real)pVert0->e[1], (real)pVert0->e[2]);
        Vec p1((real)pVert1->e[0], (real)pVert1->e[1], (real)pVert1->e[2]);
        Vec p2((real)pVert2->e[0], (real)pVert2->e[1], (real)pVert2->e[2]);
        Vec s0 = p1 - p0;
        Vec s1 = p2 - p0;
        pFaceNs[i] = (s0 % s1).normalize();
        
        for (int iVert = 0; iVert < pFace->vertex_count; iVert++) {
            int vidx = pFace->vertex_index[iVert];
            pVertNs[vidx] += pMesh->pFaces[i].normal;
        }
    }
    for (int i = 0; i < pLoader->vertexCount; i++) {
        pVertNs[i].normalize();
    }
    
    Vertex* pV = pMesh->pVertices;
    u32 iFace = 0;
    for (int i = 0; i < pLoader->faceCount; i++) {
        obj_face* pFace = pLoader->faceList[i];
        
        if (pFace->vertex_count == 3) {
            pMesh->pFaces[iFace].normal = pFaceNs[i];
            
            for ( int iVert = 0; iVert < 3; iVert++ ) {
                int vidx = pFace->vertex_index[iVert];
                pMesh->pFaces[iFace].indices[iVert] = vidx;
                
                int nidx = pFace->normal_index[iVert];
                if (nidx != -1 && nidx < pLoader->normalCount) {
                    obj_vector* pNorm = pLoader->normalList[nidx];
                    pV->normal.x = (float)pNorm->e[0];
                    pV->normal.y = (float)pNorm->e[1];
                    pV->normal.z = (float)pNorm->e[2];
                } else {
                    pV->normal = pVertNs[vidx];
                }
                
                //DebugPrint( L"%f %f %f¥n", pCv->x, pCv->y, pCv->z );
            }
            iFace++;
        }
        else if (pFace->vertex_count == 4)
        {
            // 0, 1, 2
            pMesh->pFaces[iFace].normal = pFaceNs[i];
            
            for (int iVert = 0; iVert < 3; iVert++)
            {
                int vidx = pFace->vertex_index[iVert];
                pMesh->pFaces[iFace].indices[iVert] = vidx;
                
                int nidx = pFace->normal_index[iVert];
                if (nidx != -1 && nidx < pLoader->normalCount) {
                    obj_vector* pNorm = pLoader->normalList[nidx];
                    pV->normal.x = (float)pNorm->e[0];
                    pV->normal.y = (float)pNorm->e[1];
                    pV->normal.z = (float)pNorm->e[2];
                } else {
                    pV->normal = pVertNs[vidx];
                }
                
                pV++;
            }
            iFace++;
            
            // 0, 2, 3
            int vertLocalIndices[] = { 0, 2, 3 };
            pMesh->pFaces[iFace].normal = pFaceNs[i];
            
            for (int iVert = 0; iVert < 3; iVert++) {
                int vertLocalIdx = vertLocalIndices[iVert];
                int vidx = pFace->vertex_index[vertLocalIdx];
                pMesh->pFaces[iFace].indices[iVert] = vidx;
                
                int nidx = pFace->normal_index[vertLocalIdx];
                if (nidx != -1 && nidx < pLoader->normalCount) {
                    obj_vector* pNorm = pLoader->normalList[nidx];
                    pV->normal.x = (float)pNorm->e[0];
                    pV->normal.y = (float)pNorm->e[1];
                    pV->normal.z = (float)pNorm->e[2];
                } else {
                    pV->normal = pVertNs[vidx];
                }
                
                pV++;
            }
            iFace++;
        }
        else
        {
            assert( false );
        }
    }
    
    return pMesh;
}
