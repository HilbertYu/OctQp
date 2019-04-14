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

    template <int Dim>
    class RnData
    {
        std::vector<double> pt;
    public:
        RnData(void):
            pt(Dim, 0)
        {
            assert(Dim >= 0);
        };

        const double & operator()(int idx) const
        {
            return  pt.at(idx);
        }

        double & operator()(int idx)
        {
            return  pt.at(idx);
        }

        int dim(void) const
        {
            return Dim;
        }

    };

    typedef RnData<2> R2Data;

    template <int Dim>
    class RnDataList
    {
        typedef RnData<Dim> iDataType;
        std::vector<iDataType> data_list;
    public:

        void push_back(const iDataType  & data)
        {
            data_list.push_back(data);
        }

        void clear(void)
        {
            data_list.clear();

        }

        size_t size(void) const
        {
            return data_list.size();
        }

        iDataType & operator[](int idx)
        {
            return data_list.at(idx);
        }

        const iDataType & operator[](int idx) const
        {
            return data_list.at(idx);
        }


    };

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
        RnDataList<2> data_list;

    public:

        template<typename pos_t>
        RegressionLine(const pos_t & vx, const pos_t & vy, int n_data)
        {
            data_list.clear();
            for (int i = 0; i < n_data; ++i)
            {
                R2Data v;
                v(0) = vx[i];
                v(1) = vy[i];
                data_list.push_back(v);
            }
        }

        LineEq lsm_line(void) const;

        LineEq max_norm_line(void) const;


    };
}

#endif /* end of include guard: HYOCT_H */

