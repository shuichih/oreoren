#ifndef StringUtils_H
#define StringUtils_H

#include <vector>
#include <string>

using namespace std;

class StringUtils
{
private:
    StringUtils();
    
public:
    static bool IsTrimChar(char ch);
    static string Trim(const string& str);
    static char* Trim(char* pStr);
    static int Find(const char* pStr, char ch);
    static vector<string> Split(const string& str, char delim);
    static char* Split(char* pStr, char delim);
    static int Stricmp(string lhs, string rhs);
    static int Sprintf(char* pOut, size_t szBuf, const char* pFormat, ...);
    static bool ParseFloat(float* pOut, const char* pStr);
    static bool ParseInt(int* pOut, const char* pStr);
};

#endif
