//
//  App.h
//  rt
//
//  Created by 秀一 林 on 11/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef rt_App_h
#define rt_App_h

class Photon_map;

class App
{
public:
    enum Mode
    {
        SimpleRayTrace,
        PhotonMapping
    };

    App();
    void Run(int argc, const char * argv[], int w, int h, Mode mode);
    void Init(int argc, const char * argv[], int w, int h, Mode mode);
    void Update();
    
private:
    int w_;
    int h_;
    Mode mode_;
    Photon_map* pPhotonMap_;
//    unsigned int* pColorBuf_;
};

#endif
