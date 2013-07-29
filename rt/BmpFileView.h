#ifndef Rayzer_BmpFileView_H
#define Rayzer_BmpFileView_H

#include "Common.h"
#include "IView.h"

#ifdef _WIN32

/**
 * BitMap形式でファイルに書き出すView
 */
class BmpFileView : public IView
{
public:
    
    BmpFileView();
    virtual ~BmpFileView();
    virtual bool Init(i32 width, i32 height);
    virtual bool Present(u8* pColorBuf);

private:

    bool initialized_;
    int width_;
    int height_;
};

#endif // _WIN32

#endif // Rayzer_BmpFileView_H
