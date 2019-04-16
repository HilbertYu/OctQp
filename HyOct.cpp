#include "HyOct.h"
#include <thread>

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


 //   template <typename Func_t, typename T>
    //void repeatMed(T & vec, Func_t func, const RnDataList<2> &dl)
    //void repeatMed(T & vec, Func_t func)
    
    class TSRCtxBase
    {
    public:
        virtual double cal_func(const HyOct::RnDataList<2> & dl, int i, int j) = 0;

        double repeatMed(std::vector<double> & vec, const HyOct::RnDataList<2> & dl)
        {
            const int N = dl.size();

            for (int i = 0; i < N; ++i)
            {
                std::vector<double> tmp_col_s(N, 0);
                for (int j = 0; j < N; ++j)
                {
                    tmp_col_s[j] = (cal_func(dl, i, j));
                }

                nth_element(tmp_col_s.begin(), tmp_col_s.begin() + N/2,  tmp_col_s.end());
                vec[i] = (tmp_col_s.at(N/2));
            }

            nth_element(vec.begin(), vec.begin() + N/2,  vec.end());
            return vec[N/2];

        }
    };

    double repeatMed(TSRCtxBase & b, std::vector<double> & vec, const HyOct::RnDataList<2> &dl)
    {
        return b.repeatMed(vec, dl);
    }


    LineEq RegressionLine::tsr_line(void) const
    {
        const RnDataList<2> & dl = data_list;

        auto slop = [](const RnDataList<2> & dl, int i, int j)->double
        {
            double ret = 0;
            double dem  = dl[j](0) - dl[i](0);

            if (dem == 0)
                return 1024;

            ret = (dl[j](1) - dl[i](1))/dem;
            return ret;
        };

        auto intp = [](const RnDataList<2> & dl, int i, int j)->double
        {
            double ret = 0;
            double dem  = dl[j](0) - dl[i](0);

            if (dem == 0)
                return 1024;

            ret = (dl[j](0)*dl[i](1) - dl[i](0)*dl[j](1))/dem;
            return ret;
        };


        const int N = data_list.size();

        using namespace std;
        vector<double> col_s(N, 0);
        vector<double> col_i(N, 0);


        class TSRCtxSlop: public TSRCtxBase
        {
        public:
            double cal_func(const HyOct::RnDataList<2> & dl, int i, int j)
            { 
                double ret = 0;
                double dem  = dl[j](0) - dl[i](0);

                if (dem == 0)
                    return 1024;

                ret = (dl[j](1) - dl[i](1))/dem;
                return ret;
            }
        };

        class TSRCtxIntp: public TSRCtxBase
        {
        public:
            double cal_func(const HyOct::RnDataList<2> & dl, int i, int j)
            { 
                double ret = 0;
                double dem  = dl[j](0) - dl[i](0);

                if (dem == 0)
                    return 1024;

                ret = (dl[j](0)*dl[i](1) - dl[i](0)*dl[j](1))/dem;
                return ret;
            }
        };


        TSRCtxSlop tslp;
        TSRCtxIntp titp;

        thread th_slop(repeatMed, std::ref(tslp),  std::ref(col_s), std::cref(dl));
        thread th_intp(repeatMed, std::ref(titp),  std::ref(col_i), std::cref(dl));

        th_slop.join();
        th_intp.join();


        HyOct::LineEq ret_eq = HyOct::LineEq(col_s[N/2], -1,  col_i[N/2]);
        return ret_eq;

    }

};
