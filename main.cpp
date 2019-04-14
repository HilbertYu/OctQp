#include <iostream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "HyOct.h"

//========= sample =======
typedef struct
{
    int coord[2];
} point_ctx_t;

typedef struct
{
    point_ctx_t pt[32];
    int n_pts;
} point_queue_t;
//===========================


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
    point_queue_t pq;

    double x[3] = {1, 3, 5};
    double y[3] = {1, 1, 2};

    {
        for (int i = 0; i < 3; ++i)
        {
            pq.pt[i].coord[0] = x[i];
            pq.pt[i].coord[1] = y[i];
        }
        pq.n_pts = 3;
    }

    auto init_f = [](int i, const point_queue_t & apq)
        ->HyOct::RnData<2>
    {
        HyOct::RnData<2> ret;
        ret(0) = apq.pt[i].coord[0];
        ret(1) = apq.pt[i].coord[1];
        return ret;
    };


    HyOct::RegressionLine mrl(init_f, pq, pq.n_pts);

    //AX + BY + C = 0;
    std::cout << mrl.max_norm_line();
    std::cout << mrl.lsm_line();

    using namespace std;

    HyOct::LineError err_ret =
        err_cal2(mrl.dataList(), mrl.max_norm_line());

    printf("max = %.6lf, mean = %.6lf, rms = %.6lf\n",
        err_ret.norm_max,
        err_ret.norm_mean,
        err_ret.rms);

    return 0;
}
