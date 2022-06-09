//
// Created by Flower on 2020/8/28.
//

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/owens_t.hpp>

using std::cout;
using std::endl;

double normalCDF(double x) {
    return std::erfc(-x / std::sqrt(2)) / 2;
}

double normalPDF(const double & x) {
    return (exp(-0.5 * x * x) / (double)pow(2 * M_PI, 0.5));
}

int main() {
    double a = 2;
    cout.precision(50);
    auto ori = [a] (double h) {return normalCDF(h) - 2 * boost::math::owens_t(h, a);};
    auto APS_UT = [] (double h) {return 2 * normalCDF(h) - 1;};
    auto AZC_UT = [] (double h) {return 1 - 2 * normalPDF(h) / h;};
    auto APS_LT = [a] (double h) {return 2 * normalCDF(h) * normalCDF(a * h) / (1 + a * a);};
    auto AZC_LT = [a] (double h) {return std::sqrt(2 / M_PI) * normalPDF(h * std::sqrt(1 + a * a)) / (a * (1 + a * a) * h * h);};

    double m1 = 1.0000000;
    double m2 = 1 - 1e-17;
    cout << m1 << " | " << m2 << endl;

    if (m2 == 1) {
        cout << "[1 - 1e-17] equals 1 for computers!" << endl;
    }

    return 0;
}

// g++ -std=c++11 main.cpp -o main
// g++ -bundle -DSYSTEM=APPLEMAC stplugin.c CSkewNormal.cpp -o CSkewNormal.plugin