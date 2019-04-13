#include <iostream>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "HyOct.h"

int main (void)
{

    double x[3] = {1.5, 3, 5};
    double y[3] = {1, 1, 2};
    HyOct::MaxRegressionLine mrl(x, y, 3);
    mrl.TestShow();
    std::vector<double> r = mrl.run();
    for (int i = 0; i < 3; ++i)
        printf("r(%d) = %.2lf\n", i, r[i]);

    return 0;
}
