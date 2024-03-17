//
//  main.cpp
//
//
//  Created by Richard Edberg on 2017-10-10.
//  Copyright © 2017 Richard Edberg. All rights reserved.
//

#include <unistd.h> //in order to use getcwd
#include <iostream>
#include <stdlib.h>
#include "JobStarter.h"

//========================================================================


std::string ExePath() {
    //This void returns the path where the exec file is located
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) == NULL)
    {
        std::cerr<<"Error with finding current directory";
        exit(EXIT_FAILURE);
    }
    char *getcwd(char *buf, size_t size);
    return getcwd(cwd, sizeof(cwd));
}


int main(int argc, char *argv[])
{
    bool automaticPath= false;
    std::string pathToDataFolder="";
    if(automaticPath)
    {
        pathToDataFolder= ExePath()+"/";
    }else
    {
        pathToDataFolder= "/Users/richardedberg/Documents/C++/Pyrochlore/";
    }
    
    
    JobStarter jobStarter(argc, argv, pathToDataFolder);
    
    jobStarter.startTrivialMPI();
    
    
    
    
    
    
    
    
    
}

