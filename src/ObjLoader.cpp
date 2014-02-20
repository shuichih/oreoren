#include "ObjLoader.h"
#include "File.h"
#include "StringUtils.h"
#include <cstdio>
#include <cassert>

namespace
{

enum Key
{
    Key_Invalid = -1,
    Key_V = 0,
    Key_VN,
    Key_F,
};
    
static const char* s_Keywords[] = {
    "v",
    "vn",
    "f"
};

} // namespace anonymous

bool (ObjLoader::*s_pParseFuncs[3])(char*);

ObjLoader::ObjLoader()
{
    s_pParseFuncs[Key_V]  = &ObjLoader::ParseV;
    s_pParseFuncs[Key_VN] = &ObjLoader::ParseVN;
    s_pParseFuncs[Key_F]  = &ObjLoader::ParseF;
}

ObjLoader::~ObjLoader()
{
}

bool ObjLoader::Load(const char* pFilePath)
{
    File f;
    if (!f.Open(pFilePath, "r")) {
        return false;
    }
    static char line[1024];
    int nLine = 0;
    for ( ;!f.IsEof(); nLine++) {
        line[0] = '\0';
        f.GetLine(line, 1024);
        //printf(line);
        char* pLine = StringUtils::Trim(line);
        //char* pLine = line;
        if (pLine == NULL) {
            continue;
        }
        char* pKey = pLine;
        char* pValue = StringUtils::Split(pLine, ' ');
        if (pValue == pLine) {
            continue;
        }
        pValue = StringUtils::Trim(pValue);
        int key = GetKey(pKey);
        if (key == Key_Invalid) {
            continue;
        }
        if (!(this->*s_pParseFuncs[key])(pValue)) {
            printf("Obj file parse error! line=%d\n", nLine);
        }
    }
    nVertices = (int)vertices_.size();
    nNormals = (int)normals_.size();
    nFaces = (int)faces_.size();
    return true;
}

void ObjLoader::Reset()
{
    faces_.clear();
    vertices_.clear();
    normals_.clear();
}

int ObjLoader::GetKey(char* pKeyword)
{
    for (int i=0; i<ARRAY_SZ(s_Keywords); i++) {
        if (strcmp(pKeyword, s_Keywords[i]) == 0) {
            return (Key)i;
        }
    }
    return Key_Invalid;
}

bool ObjLoader::ParseVector(Vec3& v, char* pValue)
{
    char* pX = pValue;
    char* pY = StringUtils::Trim(StringUtils::Split(pX, ' '));
    char* pZ = StringUtils::Trim(StringUtils::Split(pY, ' '));
    if (pX == pY || pY == pZ) return false;
    
    if (!StringUtils::ParseFloat(&v.x, pX)) return false;
    if (!StringUtils::ParseFloat(&v.y, pY)) return false;
    if (!StringUtils::ParseFloat(&v.z, pZ)) return false;
    
    return true;
}

bool ObjLoader::ParseVTN(int& v, int& t, int& n, char* pValue)
{
    char* pV = pValue;
    char* pT = StringUtils::Split(pV, '/');
    char* pN = StringUtils::Split(pT, '/');
    
    if (!StringUtils::ParseInt(&v, pV)) return false;
    
    t = 0;
    n = 0;
    if (pV != pT) {
        if (!StringUtils::ParseInt(&t, pT)) return false;
    } else if (pT != pN) {
        if (!StringUtils::ParseInt(&n, pN)) return false;
    }
    
    // indexが1始まりなので
    v--;
    t--;
    n--;
    
    return true;
}

bool ObjLoader::ParseV(char* pValue)
{
    Vec3 v;
    if (!ParseVector(v, pValue)) {
        return false;
    }
    vertices_.push_back(v);
    return true;
}

bool ObjLoader::ParseVN(char* pValue)
{
    Vec3 v;
    if (!ParseVector(v, pValue)) return false;
    normals_.push_back(v);
    return true;
}

bool ObjLoader::ParseF(char* pValue)
{
    char* pV0 = pValue;
    char* pV1 = StringUtils::Trim(StringUtils::Split(pV0, ' '));
    char* pV2 = StringUtils::Trim(StringUtils::Split(pV1, ' '));
    char* pV3 = StringUtils::Trim(StringUtils::Split(pV2, ' '));
    if (pV0 == pV1 || pV1 == pV2) return false;
    
    ObjFace f;
    if (!ParseVTN(f.iVertex[0], f.iTexture[0], f.iNormal[0], pV0)) return false;
    if (!ParseVTN(f.iVertex[1], f.iTexture[1], f.iNormal[1], pV1)) return false;
    if (!ParseVTN(f.iVertex[2], f.iTexture[2], f.iNormal[2], pV2)) return false;
    if (pV2 == pV3) {
        f.nVertices = 3;
    } else {
        if (!ParseVTN(f.iVertex[3], f.iTexture[3], f.iNormal[3], pV2)) return false;
        f.nVertices = 4;
    }
    
    faces_.push_back(f);
    
    return true;
}
