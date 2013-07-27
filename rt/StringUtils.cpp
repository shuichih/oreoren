#include "Common.h"
#include "StringUtils.h"
#include <algorithm>
#include <string>
#include <ctype.h>

//

string StringUtils::Trim(const string& str)
{
    size_t start = 0;
    size_t len = str.length();
    for (size_t i=0; i<str.length(); i++) {
        if (str[i] == ' ' || str[i] == '\t' || str[i] == '\n')
            start++;
        else
            break;
    }
    for (size_t i=str.length()-1; i>=0; i--) {
        if (str[i] == ' ' || str[i] == '\t' || str[i] == '\n')
            len--;
        else
            break;
    }
    len -= start;
    
    return str.substr(start, len);
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

int StringUtils::Stricmp(string lhs, string rhs)
{
    std::transform(lhs.begin(), lhs.end(), lhs.begin(), tolower);
    std::transform(rhs.begin(), rhs.end(), rhs.begin(), tolower);
    return lhs.compare(rhs);
}
