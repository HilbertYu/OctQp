#include <iostream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "HyOct.h"
#include <string>
#include <stdio.h>
#include <fstream>
#include <thread>

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


HyOct::RegressionLine TestInitRLByFunction(void)
{
    using namespace HyOct;
    point_queue_t pq;
    TestInitPQ(pq);

    //init callback function
    auto init_f = [](int i, const point_queue_t & apq)
        ->RnData<2>
    {
        RnData<2> ret;
        ret(0) = apq.pt[i].coord[0];
        ret(1) = apq.pt[i].coord[1];
        return ret;
    };

    return RegressionLine(init_f, pq, pq.n_pts);
}


void TestShowInfo(const std::vector<HyOct::LineEq> & lines, const HyOct::RnDataList<2> & rn_data)
{
    using namespace std;
    using namespace HyOct;

    vector<string> ts;
    ts.push_back("max-line");
    ts.push_back("lsm-line");
    ts.push_back("tsr-line");

    for (int i: {0,1, 2})
    {
        printf("=== %s ===\n", ts[i].c_str());

        LineError el = LineError::calError(rn_data, lines[i]);

        cout << el << endl;
    }

}


int main(int argc, const char *argv[])
{
    //Init PQ to demo


    using namespace HyOct;
    using namespace std;



    const char * file_name = "pts";
    RnDataList<2> file_data = RnDataList<2>::R2DataFileLoader(file_name);
    RegressionLine rl(file_data);

    //get the lines
    LineEq tsr_line = rl.tsr_line();
    LineEq max_line = rl.max_norm_line();
    LineEq lsm_line = rl.lsm_line();


    //(A,B,C) means AX + BY + C = 0
    cout << "max-line eq : " << max_line;
    cout << "lsm-line eq : " << lsm_line;
    cout << "tsr-line eq : " << tsr_line;
    cout << endl;

    //get the errors
    const RnDataList<2> & rn_data = rl.dataList();


    vector<LineEq> lines;

    lines.push_back(max_line);
    lines.push_back(lsm_line);
    lines.push_back(tsr_line);

    TestShowInfo(lines, rn_data);





    return 0;
}
