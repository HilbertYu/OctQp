#ifndef HYOCT_H
#define HYOCT_H

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/builtin-defun-decls.h>
#include <octave/CMatrix.h>
#include <octave/Array.h>


namespace HyOct
{

    template<typename T>
    ColumnVector array2col(const T & a, int n)
    {
        ColumnVector ret(n);
        for (int i = 0; i < n; ++i)
            ret(i) = a[i];
        return ret;
    }


    class OctQP
    {
        ColumnVector m_x0;
        Matrix m_H;
        ColumnVector m_q;
        Matrix m_Aeq;
        ColumnVector m_Beq;
        Matrix m_Ain;
        ColumnVector m_Bin;
    public:
        OctQP(const ColumnVector & x0,
           const Matrix & H,
           const ColumnVector & q,
           const Matrix & Aeq,
           const ColumnVector & Beq,
           const Matrix & Ain,
           const ColumnVector & Bin);

        ColumnVector minimize(void);
    };

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




    class RegressionLine
    {
        typedef std::vector<double> rn_data;
        std::vector<rn_data> data_list;

    public:

        template<typename pos_t>
        RegressionLine(const pos_t & vx, const pos_t & vy, int n_data)
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

        LineEq lsm_line(void) const;

        LineEq max_norm_line(void) const;


    };
}

#endif /* end of include guard: HYOCT_H */
