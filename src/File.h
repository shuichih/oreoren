#ifndef Rayzer_File_h
#define Rayzer_File_h

#include "Common.h"
#include <string>

/**
 * File access class
 */
class File
{
public:
    File();
    File(const char* pFilePath);
    ~File();
    
    bool Open(const char* pFilePath, const char* pMode);
    bool Close();
    bool IsEof();
    char* GetLine(char* pBuf, int szBuf);
    std::string GetLineS();
    std::string GetErrorString();
    
private:
    FILE* fp;
};

#endif