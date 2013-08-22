#include "BmpFileView.h"
#include "bmpexporter/bmpexporter.h"
#include <cstdio>
#include <string>
#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

using namespace std;


BmpFileView::BmpFileView()
: initialized_(false)
{
}

BmpFileView::~BmpFileView()
{
}

bool BmpFileView::Init(int w, int h)
{
    if (initialized_) {
        return true;
    }

    width_ = w;
    height_ = h;

    initialized_ = true;

    return true;
}

bool BmpFileView::Present(u8* pColorBuf)
{
    if (!initialized_) {
        printf("The view is not initialized.\n");
        return false;
    }

    string fnm = "./rendered_image.bmp";

    // Convert to 24bpp
    u8* p24bppBuf = new u8[width_*height_*3];
    for (int i = 0, j=0; i < (width_*height_*3); i+=3, j+=4)
    {
        p24bppBuf[i+0] = pColorBuf[j+0];
        p24bppBuf[i+1] = pColorBuf[j+1];
        p24bppBuf[i+2] = pColorBuf[j+2];
    }

    int result = exportToBmp(fnm.c_str(), p24bppBuf, width_, height_);
    if (result == 0) {
        printf("Failed to export bmp.\n");
        delete p24bppBuf;
        return false;
    }
    delete p24bppBuf;

#ifdef _WIN32
    ShellExecuteA(NULL, "open", fnm.c_str(), NULL, NULL, SW_SHOWNORMAL);
#endif


    return true;
}

