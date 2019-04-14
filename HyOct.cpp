#include "HyOct.h"

namespace HyOct
{

    //=============Line Error======


    std::ostream &operator<<(std::ostream &s, const LineError & le)
    {
        printf("max = %.6lf, mean = %.6lf, rms = %.6lf\n",
                le.norm_max,
                le.norm_mean,
                le.rms);
        return s;
    }

    //=========================

    OctQP::OctQP(const ColumnVector & x0,
            const Matrix & H,
            const ColumnVector & q,
            const Matrix & Aeq,
            const ColumnVector & Beq,
            const Matrix & Ain,
            const ColumnVector & Bin):
        m_x0(x0),
        m_H(H),
        m_q(q),
        m_Aeq(Aeq),
        m_Beq(Beq),
        m_Ain(Ain),
        m_Bin(Bin)
    {

    }

    ColumnVector OctQP::minimize(void)
    {
        octave_value_list x;

        x(0) = m_x0;
        x(1) = m_H;
        x(2) = m_q;
        x(3) = m_Aeq;
        x(4) = m_Beq;
        x(5) = m_Ain;
        x(6) = m_Bin;
        x(7) = 200;

        octave_value_list ret = F__qp__(x, 1);
        return ret(0).vector_value();
    }

    std::ostream &operator<<(std::ostream &s, const LineEq & line)
    {
        printf("%.3lf, %.3lf, %.3lf\n",
                line.a, line.b, line.c);
        return s;
    }

    //========= RegressionLine

    LineEq RegressionLine::lsm_line(void) const
    {
        const int n_data = data_list.size();

        Matrix A(n_data, 2);
        ColumnVector b(n_data);
        for (int i = 0; i < n_data; ++i)
        {
            A(i, 0) = data_list[i](0);
            A(i, 1) = 1;
            b(i) = data_list[i](1);
        }

        ColumnVector w = A.solve(b);
        return LineEq(w(0), -1, w(1));
    }

    LineEq RegressionLine::max_norm_line(void) const
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
                T(i, axis) = data_list[i](axis);
            }
        }

        T = T.append(ColumnVector(n_data).fill(1));
        T = T.stack(-T);
        Ain = T;

        Bin.fill(-1);

        ColumnVector r 
            = OctQP(x0, H, q, Aeq, Beq, Ain, Bin).minimize();

        return LineEq(r(0), r(1), r(2));
    }


};
