#include "LineInfo.h"

#include <stdlib.h>

namespace HyOct
{
    LineError 
    LineError::calError(const RnDataList<2> & dl, const LineEq & line_eq)
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
}

