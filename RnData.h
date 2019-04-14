#ifndef RNDATA_H
#define RNDATA_H

#include <vector>
#include <assert.h>
#include <iostream>
#include <math.h>

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
            const RnData & x = *this;
            return  const_cast<double&>(x(idx));
        }

        int dim(void) const
        {
            return Dim;
        }

    };


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

};

#endif /* end of include guard: RNDATA_H */

