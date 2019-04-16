#include <iostream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "HyOct.h"
#include <string>
#include <stdio.h>
#include <fstream>

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

HyOct::RnDataList<2> R2DataFileLoader(const std::string & file_name)
{
    using namespace HyOct;
    using namespace std;

    RnDataList<2> ret;

    ifstream ifs(file_name);
    if (ifs.bad())
    {
        fprintf(stderr, "error\n");
        exit(-1);
    }

    string line;
    while (getline(ifs, line))
    {
        if (ifs.eof())
            break;

        if (line.size() == 0)
            break;

        int coord[2] = {-1, -1};
        int r = sscanf(line.c_str(), "%d,%d\n", coord, coord+1);
        assert(r == 2);

        RnData<2> v;
        v(0) = coord[0];
        v(1) = coord[1];
        ret.push_back(v);
    }

    return ret;

}



int main(int argc, const char *argv[])
{
    //Init PQ to demo
    point_queue_t pq;
    TestInitPQ(pq);


    using namespace HyOct;
    using namespace std;

    //init callback function
    auto init_f = [](int i, const point_queue_t & apq)
        ->RnData<2>
    {
        RnData<2> ret;
        ret(0) = apq.pt[i].coord[0];
        ret(1) = apq.pt[i].coord[1];
        return ret;
    };


    const char * file_name = "pts";
    RnDataList<2> file_data = R2DataFileLoader(file_name);
    {
        for (size_t i = 0; i < file_data.size(); ++i)
        {
            double x = file_data[i](0);
            double y = file_data[i](1);
            printf("PT,%.3lf, %.3lf\n", x, y);

        }
    }

    //RegressionLine rl(init_f, pq, pq.n_pts);
    RegressionLine rl(file_data);

    //get the lines
    LineEq max_line = rl.max_norm_line();
    LineEq lsm_line = rl.lsm_line();


    //(A,B,C) means AX + BY + C = 0
    cout << "max-line eq : " << max_line;
    cout << "lsm-line eq : " << lsm_line;

    //get the errors
    const RnDataList<2> & rn_data = rl.dataList();

    LineError err_ret_max_line =
        LineError::calError(rn_data, max_line);

    LineError err_ret_lsm_line =
        LineError::calError(rn_data, lsm_line);

    cout << err_ret_max_line;
    cout << err_ret_lsm_line;

    return 0;
}
