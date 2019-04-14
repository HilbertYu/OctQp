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

HyOct::LineError
err_cal2(
        const HyOct::RnDataList<2> & dl,
        const HyOct::LineEq & line_eq)
{

    const double A = line_eq.a;
    const double B = line_eq.b;
    const double C = line_eq.c;

    int n_data = dl.size();

    double w = sqrt(A*A + B*B);

    double ret_max_err = -1;
    double ret_mean_err = 0;
    double ret_rms_err = 0;


    for (int i = 0; i < n_data; ++i)
    {
        const HyOct::RnData<2> & pt = dl[i];
        double r = fabs(A*pt(0) + B*pt(1) + C)/w;
        if (ret_max_err < r)
            ret_max_err = r;

        ret_mean_err += r;
        ret_rms_err += r*r;

    }

    HyOct::LineError ret;
    ret.norm_max = ret_max_err;
    ret.norm_mean = ret_mean_err/n_data;
    ret.rms = sqrt(ret_rms_err/n_data);

    return ret;

}

int main(int argc, const char *argv[])
{

    double x[3] = {1, 3, 5};
    double y[3] = {1, 1, 2};

  //  RegressionLine(InitFunc func, const ListT & dl, int n_data);

    class EP
    {
    public:
        double * x;
        double * y;
    };

    EP ep;
    ep.x = x;
    ep.y = y;

    auto init_f = [](int i, const EP & p)->HyOct::RnData<2>
    {
        HyOct::RnData<2> ret;
        ret(0) = p.x[i];
        ret(1) = p.y[i];
        return ret;
    };


//    HyOct::RegressionLine mrl(x, y, 3);
    HyOct::RegressionLine mrl(init_f, ep, 3);

    //AX + BY + C = 0;

    std::cout << mrl.max_norm_line();
    std::cout << mrl.lsm_line();

    using namespace std;
    cout << err_cal(x, y, 3, mrl.max_norm_line()) << endl;

    HyOct::LineError err_ret =
        err_cal2(mrl.dataList(), mrl.max_norm_line());

    printf("max = %.6lf, mean = %.6lf, rms = %.6lf\n",
        err_ret.norm_max,
        err_ret.norm_mean,
        err_ret.rms);

    return 0;
}
