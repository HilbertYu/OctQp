#ifndef HYOCT_H
#define HYOCT_H

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/builtin-defun-decls.h>
#include <octave/CMatrix.h>
#include <octave/Array.h>

#include "RnData.h"
#include "LineInfo.h"

namespace HyOct
{
    template <int Dim> class RnData;
    template <int Dim> class RnDataList;

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
        RnDataList<2> data_list;

    public:
        template <typename InitFunc, typename ListT>
        RegressionLine(InitFunc func, const ListT & dl, int n_data)
        {
            data_list.clear();
            for (int i = 0; i < n_data; ++i)
            {
                RnData<2> v = func(i, dl);
                data_list.push_back(v);
            }

        }

        const RnDataList<2> & dataList(void) const
        {
            return data_list;
        }

        template <int Dim>
        RegressionLine(const RnDataList<Dim> & arg_data_list):
            data_list(arg_data_list)
        {}


        LineEq lsm_line(void) const;

        LineEq max_norm_line(void) const;

        LineEq tsr_line(void) const
        {
            const RnDataList<2> & dl = data_list;

            auto slop = [&dl](int i, int j)->double
            {
                double ret = 0;
                double dem  = dl[j](0) - dl[i](0);

                if (dem == 0)
                    return 1024;

                ret = (dl[j](1) - dl[i](1))/dem;
                return ret;
            };

            auto intp = [&dl](int i, int j)->double
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


            for (int i = 0; i < N; ++i)
            {
                vector<double> tmp_col_s(N, 0);
                vector<double> tmp_col_i(N, 0);

                for (int j = 0; j < N; ++j)
                {

                    tmp_col_s[j] = (slop(i, j));
                    tmp_col_i[j] = (intp(i, j));
                }

                nth_element(tmp_col_s.begin(), tmp_col_s.begin() + N/2,  tmp_col_s.end());
                nth_element(tmp_col_i.begin(), tmp_col_i.begin() + N/2,  tmp_col_i.end());
                
                col_s[i] = (tmp_col_s.at(N/2));
                col_i[i] = (tmp_col_i.at(N/2));
            }


            // for (int i = 0; i < N; ++i)
            // {
            //     printf("T,%.3lf\n", col_i[i]);
            //
            // }

            nth_element(col_s.begin(), col_s.begin() + N/2,  col_s.end());
            nth_element(col_i.begin(), col_i.begin() + N/2,  col_i.end());


            HyOct::LineEq ret_eq = HyOct::LineEq(col_s[N/2], -1,  col_i[N/2]);
            return ret_eq;

        }


    };
}

#endif /* end of include guard: HYOCT_H */

