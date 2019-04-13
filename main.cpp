#include <iostream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "HyOct.h"

int main(int argc, const char *argv[])
{

    double x[3] = {1, 3, 5};
    double y[3] = {1, 1, 2};
    HyOct::RegressionLine mrl(x, y, 3);

    std::vector<double> r = mrl.max_norm_line();
    r = mrl.lsm_line();

    for (int i = 0; i < 3; ++i)
        printf("r(%d) = %.2lf\n", i, r[i]);


    Matrix A(3,2);
    A.fill(1);
    ColumnVector b(3);

    for (int i = 0; i < 3; ++i)
    {
        A(i, 0) = x[i];
        b(i) = y[i];
    }


    ColumnVector v = A.solve(b);
    std::cout << A;
    std::cout << b;
    std::cout << v;

    return 0;
}
