//
// Created by Flower on 2020/8/28.
//

#include "stplugin.h"

STDLL stata_call(int argc, char *argv[])
{
    SF_display("Hello World\n") ;
    return(0) ;
}