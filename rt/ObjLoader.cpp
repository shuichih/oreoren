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
    char line[1024];
    int nLine = 0;
    while (!f.IsEof()) {
        f.GetLine(line, 1024);
        //printf(line);
        char* pLine = StringUtils::Trim(line);
        if (pLine == NULL) {
            continue;
        }
        int iSplit = StringUtils::Find(pLine, ' ');
        if (iSplit == -1) {
            continue;
        }
        pLine[iSplit] = NULL;
        char* pValue = StringUtils::Trim(&pLine[iSplit+1]);
        int key = GetKey(pLine);
        if (key == Key_Invalid) {
            continue;
        }
        if (!(this->*s_pParseFuncs[key])(pValue)) {
            printf("Obj file parse erroro! line=%d\n", nLine);
        }
        nLine++;
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
    char* pY = StringUtils::Split(pX, ' ');
    char* pZ = StringUtils::Split(pY, ' ');
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
    if (pV == pT || pT == pN) return false;
    
    if (!StringUtils::ParseInt(&v, pV)) return false;
    if (!StringUtils::ParseInt(&t, pT)) return false;
    if (!StringUtils::ParseInt(&n, pN)) return false;
    
    // indexが1始まりなので
    v--;
    t--;
    n--;
    
    return true;
}

bool ObjLoader::ParseV(char* pValue)
{
    Vec3 v;
    if (!ParseVector(v, pValue)) return false;
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
    char* pV1 = StringUtils::Split(pV0, ' ');
    char* pV2 = StringUtils::Split(pV1, ' ');
    char* pV3 = StringUtils::Split(pV2, ' ');
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
