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

        std::vector<double> lsm_line(void) const
        {
            const int n_data = data_list.size();

            Matrix A(n_data, 2);
            ColumnVector b(n_data);
            for (int i = 0; i < n_data; ++i)
            {
                A(i, 0) = data_list[i][0];
                A(i, 1) = 1;
                b(i) = data_list[i][1];
            }

            ColumnVector w = A.solve(b);

            std::vector<double> ret;
            ret.push_back(w(0));
            ret.push_back(-1);
            ret.push_back(w(1));

            return ret;
        }

        std::vector<double> max_norm_line(void) const
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

#endif /* end of include guard: HYOCT_H */
