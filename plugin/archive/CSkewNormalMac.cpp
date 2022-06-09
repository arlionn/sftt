//
// Created by Flower on 2020/9/5.
//

#include "stplugin.h"
#include <cmath>
#include <boost/math/special_functions/owens_t.hpp>

auto normalCDF = [] (const double & x) {return std::erfc(-x / std::sqrt(2)) / 2;};
auto normalPDF = [] (const double & x) {return (exp(-0.5 * x * x) / (double)pow(2 * M_PI, 0.5));};
auto APS_UT = [] (const double & Q) {return 2 * normalCDF(Q) - 1;};
auto AZC_UT = [] (const double & Q) {return 1 - 2 * normalPDF(Q) / Q;};
auto APS_LT = [] (const double & Q, const double & lambda) {
    return 2 * normalCDF(Q) * normalCDF(lambda * Q) / (1 + lambda * lambda);
};
auto AZC_LT = [] (const double & Q, const double & lambda) {
    return std::sqrt(2 / M_PI) * normalPDF(Q * std::sqrt(1 + lambda * lambda)) /
            (lambda * (1 + lambda * lambda) * Q * Q);
};

STDLL stata_call(int argc, char* argv[]) {
    ST_int        j, k ;
    ST_double     h, a, G;
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
            G = normalCDF(h) - boost::math::owens_t(h, a);
            // for upper-tail
            if (G > 1 - 1e-40) {
                G = APS_UT(h);
            }
            if (G < 1e-40) {
                G = APS_LT(h, a);
            }
            rc = SF_vstore(3, j, G);
            if (rc) return(rc);
        }
    }
    return(0) ;
}

