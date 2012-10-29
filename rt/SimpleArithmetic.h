#ifndef SimpleArithmetic_H
#define SimpleArithmetic_H

#include <string>

/**
 * 文字列で与えられた簡単な四則演算を含む式を計算する
 */
class SimpleArithmetic
{
private:
    SimpleArithmetic();
public:
    static bool Eval(float* pVal, const std::string& expression);
};

#endif
