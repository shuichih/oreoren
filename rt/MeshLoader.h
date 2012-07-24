#ifndef MeshLoader_h
#define MeshLoader_h

class Mesh;

class ObjLoader
{
public:
    ObjLoader();
    ~ObjLoader();
    
    Mesh* Load(const char* pFilePath);
};

#endif
