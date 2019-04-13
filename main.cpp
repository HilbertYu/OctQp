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

    //AX + BY + C = 0;

    std::cout << mrl.max_norm_line();
    std::cout << mrl.lsm_line();
 //   r = mrl.lsm_line();




    return 0;
}
