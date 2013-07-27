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
    static string Trim(const string& str);
    static vector<string> Split(const string& str, char delim);
    static int Stricmp(string lhs, string rhs);
    static int Sprintf(char* pOut, size_t szBuf, const char* pFormat, ...);
};

#endif
