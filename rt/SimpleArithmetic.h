#ifndef SimpleArithmetic_H
#define SimpleArithmetic_H

#include <string>

class SimpleArithmetic
{
private:
    SimpleArithmetic();
public:
    static bool Eval(float* pVal, const std::string& expression);
};

#endif
