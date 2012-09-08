//
//  App.h
//  rt
//
//  Created by 秀一 林 on 11/29/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#ifndef rt_App_h
#define rt_App_h

#include "Config.h"

class Photon_map;

class App
{
public:

    App();
    ~App();
    void Run(int argc, const char * argv[]);
    void Init(int argc, const char * argv[]);
    void Update();
    
private:
    
    Config config;
    Photon_map* pPhotonMap_;
    PhotonMapRenderer renderer_;
};

#endif
