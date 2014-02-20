#ifndef Rayzer_MeshLoader_h
#define Rayzer_MeshLoader_h

#include "Common.h"
#include <vector>

#define MAX_VERTEX_COUNT 4  // quad or triangle

typedef Vec3 ObjVertex;

struct ObjFace
{
	int iVertex[MAX_VERTEX_COUNT];
	int iNormal[MAX_VERTEX_COUNT];
	int iTexture[MAX_VERTEX_COUNT];
	int nVertices;
};

class ObjLoader
{
public:
    ObjLoader();
    ~ObjLoader();
    
    bool Load(const char* pFilePath);
    inline const std::vector<Vec3>& Vertices() { return vertices_; }
    inline const std::vector<Vec3>& Normals() { return normals_; }
    inline const std::vector<ObjFace>& Faces() { return faces_; }
    inline int NumVertices() { return nVertices; }
    inline int NumNormals() { return nNormals; }
    inline int NumFaces() { return nFaces; };
    
private:
    bool (ObjLoader::*s_pParseFuncs[3])(char*);
    
    void Reset();
    int GetKey(char* pKeyword);
    bool ParseVector(Vec3& v, char* pValue);
    bool ParseVTN(int& v, int& t, int& n, char* pValue);
    bool ParseV(char* pValue);
    bool ParseVN(char* pValue);
    bool ParseF(char* pValue);
    
    std::vector<Vec3> vertices_;
    std::vector<Vec3> normals_;
    std::vector<ObjFace> faces_;
    int nVertices;
    int nNormals;
    int nFaces;
};

#endif

