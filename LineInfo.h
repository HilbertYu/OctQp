#ifndef LINEINFO_H
#define LINEINFO_H

#include "RnData.h"

namespace HyOct
{
    class LineEq;

    class LineError
    {
    public:
        LineError(void):
            norm_max(-1),
            norm_mean(-1),
            rms(-1)
        {
        }

        double norm_max;
        double norm_mean;
        double rms;

        static 
        LineError calError(const RnDataList<2> & dl, const LineEq & line_eq);

    };

    std::ostream &operator<<(std::ostream &s, const LineError & le);

    class LineEq
    {
        // LineEq  a*X + b*Y + c = 0
    public:
        double a;
        double b;
        double c;

        template <typename T>
        LineEq(const T & v):
            a(v[0]),
            b(v[1]),
            c(v[2])
        {

        }

        LineEq(double _a, double _b, double _c):
            a(_a),
            b(_b),
            c(_c)
        {}


        LineEq(void):
            a(0),
            b(0),
            c(0)
        {}


       friend std::ostream &operator<<(std::ostream &s, const LineEq &);
    };

};


#endif /* end of include guard: LINEINFO_H */
