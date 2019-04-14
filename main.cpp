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
        HyOct::LineError::calError(mrl.dataList(), mrl.max_norm_line());

    std::cout << err_ret;

    return 0;
}
