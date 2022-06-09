//
// Created by Flower on 2020/8/28.
//

#include "stplugin.h"
#include <boost/math/special_functions/owens_t.hpp>
//#include "owen.h"

STDLL stata_call(int argc, char* argv[]) {
    ST_int        j, k ;
    ST_double     h, a, t, sum, z ;
    ST_retcode	  rc ;

    if(SF_nvars() != 3) {
        return(102) ;  	    /* not enough variables specified */
    }

    for (j = SF_in1(); j <= SF_in2(); j++) {
        if (SF_ifobs(j)) {
            rc = SF_vdata(1, j, &h);
            if(rc) return(rc);
            rc = SF_vdata(2, j, &a);
            if(rc) return(rc);
            t = boost::math::owens_t(h, a);
            rc = SF_vstore(3, j, t);
            if (rc) return(rc);
//            char str[60];
//            sprintf(str, "j: %d ｜h: %.4f | a: %.4f | t: %.8f" , j, h, a, t);
//            SF_display(str);
//            SF_display("\n");
        }
    }
    return(0) ;
}

