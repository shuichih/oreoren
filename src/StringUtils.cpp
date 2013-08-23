#include "Common.h"
#include "StringUtils.h"
#include <algorithm>
#include <string>
#include <ctype.h>
#include <cstdarg>
#include <cstdlib>
#include <cassert>

//

bool StringUtils::IsTrimChar(char ch)
{
    return (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r');
}

string StringUtils::Trim(const string& str)
{
    size_t start = 0;
    size_t len = str.length();
    const i32 clen = (i32)str.length();
    for (i32 i=0; i<clen; i++) {
        if (IsTrimChar(str[i]))
            start++;
        else
            break;
    }
    for (i32 i=clen-1; i>=0; i--) {
        if (IsTrimChar(str[i]))
            len--;
        else
            break;
    }
    len -= start;
    
    return str.substr(start, len);
}

char* StringUtils::Trim(char* pStr)
{
    size_t start = 0;
    size_t len = strlen(pStr);
    const i32 clen = (i32)len;
    for (i32 i=0; i<clen; i++) {
        if (IsTrimChar(pStr[i]))
            start++;
        else
            break;
    }
    for (i32 i=clen-1; i>=0; i--) {
        if (IsTrimChar(pStr[i]))
            len--;
        else
            break;
    }
    //len -= start;
    pStr[len] = NULL;
    
    assert(start >= 0 && start < 1024);
    return &pStr[start];
}
            
int StringUtils::Find(const char* pStr, char ch)
{
    for (int i=0; *pStr!='\0'; pStr++, i++) {
        if (*pStr == ch) {
            return i;
        }
    }
    return -1;
}

vector<string> StringUtils::Split(const string& str, char delim)
{
    vector<string> ret;
    
    string::size_type start = 0;
    string::size_type pos;
    do {
        pos = str.find(delim, start);
        if (pos == string::npos) {
            pos = str.length();
        }
        
        string ss = str.substr(start, pos - start);
        if (ss.length() != 0) {
            ret.push_back(ss);
        }
        
        start = pos + 1;
    
    } while(pos < str.length());
    
    return ret;
}

char* StringUtils::Split(char* pStr, char delim)
{
    int index = Find(pStr, delim);
    if (index == -1) {
        return pStr; // delimが見つからなかったら渡された文字列をそのまま返す
    }
    pStr[index] = NULL;
    return &pStr[index+1];
}

int StringUtils::Stricmp(string lhs, string rhs)
{
    std::transform(lhs.begin(), lhs.end(), lhs.begin(), tolower);
    std::transform(rhs.begin(), rhs.end(), rhs.begin(), tolower);
    return lhs.compare(rhs);
}

int StringUtils::Sprintf(char* pOut, size_t szBuf, const char* pFormat, ...)
{
    va_list args;
    int ret;
    
    va_start(args, pFormat);
#ifdef _WIN32
    ret = vsnprintf_s(pOut, szBuf, _TRUNCATE, pFormat, args);
#else
    ret = vsnprintf(pOut, szBuf, pFormat, args);
#endif
    va_end(args);
    return ret;
}

bool StringUtils::ParseFloat(float* pOut, const char* pStr)
{
    char* e = NULL;
    *pOut = (float)strtod(pStr, &e);
    if (*e != '\0') {
        return false;
    }
    return true;
}

bool StringUtils::ParseInt(int* pOut, const char* pStr)
{
    char* e = NULL;
    *pOut = (int)strtol(pStr, &e, 10);
    if (*e != '\0') {
        return false;
    }
    return true;
}

