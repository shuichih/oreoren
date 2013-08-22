#include <string>
#include "SimpleArithmetic.h"

using namespace std;

SimpleArithmetic::SimpleArithmetic()
{
}

bool SimpleArithmetic::Eval(float* pVal, const string& expression)
{
    string::size_type idx;
    
    idx = expression.find("+", 0);
    if (idx != string::npos) {
        float v0 = 0, v1 = 0;
        if (idx > 0) {
            if (!Eval(&v0, expression.substr(0, idx))) return false;
        }
        if (!Eval(&v1, expression.substr(idx+1))) return false;
        *pVal = v0 + v1;
        return true;
    }
    idx = expression.find("-", 0);
    if (idx != string::npos) {
        float v0 = 0, v1 = 0;
        if (idx > 0) {
            if (!Eval(&v0, expression.substr(0, idx))) return false;
        }
        if (!Eval(&v1, expression.substr(idx+1))) return false;
        *pVal = v0 - v1;
        return true;
    }
    idx = expression.find("/", 0);
    if (idx != string::npos) {
        float v0 = 0, v1 = 0;
        if (idx == 0) return false;
        if (!Eval(&v0, expression.substr(0, idx))) return false;
        if (!Eval(&v1, expression.substr(idx+1))) return false;
        *pVal = v0 / v1;
        return true;
    }
    idx = expression.find("*", 0);
    if (idx != string::npos) {
        float v0 = 0, v1 = 0;
        if (idx == 0) return false;
        if (!Eval(&v0, expression.substr(0, idx))) return false;
        if (!Eval(&v1, expression.substr(idx+1))) return false;
        *pVal = v0 * v1;
        return true;
    }
    
    char* e = NULL;
    *pVal = (float)strtod(expression.c_str(), &e);
    if (*e != '\0') {
        printf("arithmetic eval error. %s\n", expression.c_str());
        return false;
    }
    return true;
}