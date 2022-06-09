//
// Created by Flower on 2020/8/28.
//

#include "stplugin.h"
//#include "owen.hpp"

STDLL stata_call(int argc, char* argv[]) {

    ST_int        j, k ;
    ST_double     sum, z ;
    ST_retcode	  rc ;

    if(SF_nvars() < 2) {
        return(102) ;  	    /* not enough variables specified */
    }

    for (j = SF_in1(); j <= SF_in2(); j++) {
        if (SF_ifobs(j)) {
            sum = 0.0 ;
            for(k=1; k < SF_nvars(); k++) {
                rc = SF_vdata(k, j, &z);
                if(rc) return(rc) ;
                sum += z ;
            }
            rc = SF_vstore(SF_nvars(), j, sum);
            if (rc) return(rc) ;
        }
    }
    return(0) ;
}

