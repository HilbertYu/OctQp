#include "HyOct.h"

namespace HyOct
{
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
};
