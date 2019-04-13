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

        ColumnVector minimize(void)
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



#if 0
    x(0) = x0;
    x(1) = H;
    x(2) = q;
    x(3) = Aeq;
    x(4) = beq;
    x(5) = Ain;
    x(6) = Bin;
    x(7) = maxit;
#endif


    };
}

#endif /* end of include guard: HYOCT_H */
