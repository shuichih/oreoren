#include "File.h"
#include <cstdio>
#include <cstring>
#include <cerrno>

File::File()
: fp(NULL)
{
}

File::File(const char* pFilePath)
{
    
}

File::~File()
{
    Close();
}

bool File::Close()
{
    int ret = true;
    if (fp != NULL) {
        ret = fclose(fp);
        fp = NULL;
    }
    return ret;
}

bool File::Open(const char* pFilePath, const char* pMode)
{
    Close();
    
#ifdef _WIN64
    fopen_s(&fp, pPath, pMode);
#else
    fp = fopen(pFilePath, pMode);
#endif
    
    return (fp != NULL);
}

bool File::IsEof()
{
    if (fp == NULL) {
        return true;
    }
    return feof(fp) != 0;
}

char* File::GetLine(char* pBuf, int szBuf)
{
    if (szBuf == 0) {
        return pBuf;
    }
    if (szBuf == 1 || IsEof()) {
        pBuf[0] = NULL;
        return pBuf;
    }
    return fgets(pBuf, szBuf, fp);
}

std::string File::GetLineS()
{
    char buf[1024];
    buf[0] = 0;
    GetLine(buf, 1024);
    if (buf[0] == 0) {
        return std::string();
    }
    return buf;
}

std::string File::GetErrorString()
{
    char errStr[512];
#ifdef _WIN64
    strerror_s(errStr, 512, errno);
#else
    return strerror(errno);
#endif
    return errStr;
}

