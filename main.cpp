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


void TestInitPQ(point_queue_t & pq)
{
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

}


//===========================


int main(int argc, const char *argv[])
{
    //Init PQ to demo
    point_queue_t pq;
    TestInitPQ(pq);


    using namespace HyOct;
    using namespace std;

    auto init_f = [](int i, const point_queue_t & apq)
        ->RnData<2>
    {
        RnData<2> ret;
        ret(0) = apq.pt[i].coord[0];
        ret(1) = apq.pt[i].coord[1];
        return ret;
    };


    RegressionLine rl(init_f, pq, pq.n_pts);

    //AX + BY + C = 0;
    cout << rl.max_norm_line();
    cout << rl.lsm_line();


    LineError err_ret =
        LineError::calError(rl.dataList(), rl.max_norm_line());

    cout << err_ret;

    return 0;
}
