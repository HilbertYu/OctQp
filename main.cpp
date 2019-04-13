#include <iostream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "HyOct.h"


namespace HyOct
{
    class MaxRegressionLine
    {
        typedef std::vector<double> rn_data;
        std::vector<rn_data> data_list;

    public:
        template<typename pos_t>
        MaxRegressionLine(const pos_t & vx, const pos_t & vy, int n_data)
        {
            data_list.clear();
            for (int i = 0; i < n_data; ++i)
            {
                rn_data v;
                v.push_back(vx[i]);
                v.push_back(vy[i]);
                data_list.push_back(v);
            }
        }

        void TestShow(void)
        {
            auto f = [](rn_data d)
            {
                printf("(%.2lf, %.2lf)\n", d[0], d[1]);
            };

            for_each(data_list.begin(), data_list.end(), f);

        }


        std::vector<double> run(void) const
        {
            const int n_dim = 3;
            const int n_data = data_list.size();

            ColumnVector x0(n_dim);
            Matrix H(n_dim, n_dim);
            ColumnVector q(n_dim);
            Matrix Aeq(n_dim, n_dim);
            ColumnVector Beq(n_dim);
            Matrix Ain(n_data*2, n_dim);
            ColumnVector Bin(n_data*2);

            x0.fill(0);

            H.fill(0);
            for (int i = 0; i < n_dim - 1; ++i)
                H(i, i) = -1;

            q.fill(0);
            Aeq.fill(0);
            Beq.fill(0);

            Matrix T(n_data, 2);
            for (int i = 0; i < n_data; ++i)
            {
                for (int axis = 0; axis < n_dim - 1; ++axis)
                {
                    T(i, axis) = data_list.at(i).at(axis);
                }
            }

            T = T.append(ColumnVector(n_data).fill(1));
            T = T.stack(-T);
            Ain = T;

            Bin.fill(-1);

            ColumnVector r = HyOct::OctQP(x0, H, q, Aeq, Beq, Ain, Bin).minimize();

            std::vector<double> ret;
            ret.push_back(r(0));
            ret.push_back(r(1));
            ret.push_back(r(2));

            return ret;
        }


    };
}


void TestQP(void)
{

    octave_value_list x;

    ColumnVector x0(3);
    Matrix H(3, 3);
    const ColumnVector q(3);

    Matrix Aeq(3, 3);
    ColumnVector beq(3);

    Matrix Ain(6, 3);
    ColumnVector Bin(6);
    const int maxit = 200;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            H(i, j) = 0;

    H(0, 0) = H(1, 1) = -1;
    x0(0) = x0(1) = x0(2) = 0;

    Ain(0, 0) = 1; Ain(0, 1) = 1; Ain(0, 2) = 1;
    Ain(1, 0) = 3; Ain(1, 1) = 1; Ain(1, 2) = 1;
    Ain(2, 0) = 5; Ain(2, 1) = 2; Ain(2, 2) = 1;
    Ain(3, 0) = -1; Ain(3, 1) = -1; Ain(3, 2) = -1;
    Ain(4, 0) = -3; Ain(4, 1) = -1; Ain(4, 2) = -1;
    Ain(5, 0) = -5; Ain(5, 1) = -2; Ain(5, 2) = -1;

    for (int i = 0; i < 6; ++i)
        Bin(i) = -1;

    x(0) = x0;
    x(1) = H;
    x(2) = q;
    x(3) = Aeq;
    x(4) = beq;
    x(5) = Ain;
    x(6) = Bin;
    x(7) = maxit;

    octave_value_list ret = F__qp__(x, 3);

    ColumnVector v = ret(0).vector_value();
    std::cout << ret(0).vector_value();

    ColumnVector r = HyOct::OctQP(x0, H, q, Aeq, beq, Ain, Bin).minimize();
    std::cout << r;

}

void Test1(void)
{
    Matrix ma(2,2);
    Matrix mb(2,2);

    ma(0,0) = 1;
    ma(0,1) = 2;
    ma(1,0) = 3;
    ma(1,1) = 4;

    mb(0,0) = 5;
    mb(0,1) = 6;
    mb(1,0) = 7;
    mb(1,1) = 8;


    using namespace std;
    cout << "ma = " << endl;
    cout << ma;
    cout << "mb = " << endl;
    cout << mb;

    cout << "ma*mb = " << endl;
    cout << ma * mb;

    cout << "ma + mb = " << endl;
    cout << ma + mb;

    cout << "transpose of ma = " << endl;
    cout << ma.transpose() << endl;

    cout << ma.stack(mb) << endl;

    ColumnVector v(3);
    v(0) = 1;
    v(1) = 2;
    v(2) = 3;

    cout << "v = " << endl;
    cout << v;

    octave_value_list input;
    input(0) = v;

    octave_value_list output = Fnorm(input, 1);

    cout << "norm of v = " << output(0).double_value() << endl;
    cout << "norm of v by C = " << sqrt(1*1 + 2*2 + 3*3) << endl;



}

int main (void)
{
//    TestQP();

    double x[3] = {1.5, 3, 5};
    double y[3] = {1, 1, 2};
    HyOct::MaxRegressionLine mrl(x, y, 3);
    mrl.TestShow();
    std::vector<double> r = mrl.run();
    for (int i = 0; i < 3; ++i)
        printf("r(%d) = %.2lf\n", i, r[i]);


//    Test1();

    return 0;
}
