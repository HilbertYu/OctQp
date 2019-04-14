#include <iostream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "HyOct.h"

template <typename T>
double err_cal(const T & vx, const T & vy, int n_data, const HyOct::LineEq & line_eq)
{
    double ret_err = -1;

    const double A = line_eq.a;
    const double B = line_eq.b;
    const double C = line_eq.c;

    double w = sqrt(A*A + B*B);

    for (int i = 0; i < n_data; ++i)
    {
        double r = fabs(A*vx[i] + B*vy[i] + C)/w;
        if (ret_err < r)
            ret_err = r;
    }

    return ret_err;

}

int main(int argc, const char *argv[])
{

    double x[3] = {1, 3, 5};
    double y[3] = {1, 1, 2};


    HyOct::RegressionLine mrl(x, y, 3);

    //AX + BY + C = 0;

    std::cout << mrl.max_norm_line();
    std::cout << mrl.lsm_line();

    using namespace std;
    cout << err_cal(x, y, 3, mrl.max_norm_line()) << endl;

 //   r = mrl.lsm_line();




    return 0;
}
