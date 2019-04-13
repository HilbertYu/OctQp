#include <iostream>
#include <vector>
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/builtin-defun-decls.h>
#include <octave/CMatrix.h>
#include <octave/Array.h>


//==============

#if 0
    ColumnVector x0(3);
    Matrix H(3, 3);
    const ColumnVector q(3);

    Matrix Aeq(3, 3);
    ColumnVector beq(3);

    Matrix Ain(6, 3);
    ColumnVector Bin(6);
    const int maxit = 200;
#endif

namespace HyOct
{
    template <int n_dim, int n_eqdim, int n_indim>
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


void TestQP(void)
{

    octave_value_list x;

    ColumnVector x0(3);
    Matrix H(3, 3);
    const ColumnVector q(3);

    Matrix Aeq(3, 3);
    ColumnVector beq(3);

    Matrix Ain(6, 3);
    ColumnVector Bin(6);
    const int maxit = 200;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            H(i, j) = 0;

    H(0, 0) = H(1, 1) = -1;
    x0(0) = x0(1) = x0(2) = 0;

    Ain(0, 0) = 1; Ain(0, 1) = 1; Ain(0, 2) = 1;
    Ain(1, 0) = 3; Ain(1, 1) = 1; Ain(1, 2) = 1;
    Ain(2, 0) = 5; Ain(2, 1) = 2; Ain(2, 2) = 1;
    Ain(3, 0) = -1; Ain(3, 1) = -1; Ain(3, 2) = -1;
    Ain(4, 0) = -3; Ain(4, 1) = -1; Ain(4, 2) = -1;
    Ain(5, 0) = -5; Ain(5, 1) = -2; Ain(5, 2) = -1;

    for (int i = 0; i < 6; ++i)
        Bin(i) = -1;

    x(0) = x0;
    x(1) = H;
    x(2) = q;
    x(3) = Aeq;
    x(4) = beq;
    x(5) = Ain;
    x(6) = Bin;
    x(7) = maxit;

    octave_value_list ret = F__qp__(x, 3);

    ColumnVector v = ret(0).vector_value();
    std::cout << ret(0).vector_value();

    ColumnVector r = HyOct::OctQP<1,1,1>(x0, H, q, Aeq, beq, Ain, Bin).minimize();
    std::cout << r;

}

void Test1(void)
{
    Matrix ma(2,2);
    Matrix mb(2,2);

    ma(0,0) = 1;
    ma(0,1) = 2;
    ma(1,0) = 3;
    ma(1,1) = 4;

    mb(0,0) = 5;
    mb(0,1) = 6;
    mb(1,0) = 7;
    mb(1,1) = 8;


    using namespace std;
    cout << "ma = " << endl;
    cout << ma;
    cout << "mb = " << endl;
    cout << mb;

    cout << "ma*mb = " << endl;
    cout << ma * mb;

    cout << "ma + mb = " << endl;
    cout << ma + mb;

    cout << "transpose of ma = " << endl;
    cout << ma.transpose() << endl;

    cout << ma.stack(mb) << endl;

    ColumnVector v(3);
    v(0) = 1;
    v(1) = 2;
    v(2) = 3;

    cout << "v = " << endl;
    cout << v;

    octave_value_list input;
    input(0) = v;

    octave_value_list output = Fnorm(input, 1);

    cout << "norm of v = " << output(0).double_value() << endl;
    cout << "norm of v by C = " << sqrt(1*1 + 2*2 + 3*3) << endl;



}

int main (void)
{
    TestQP();

//    Test1();

    return 0;
}
