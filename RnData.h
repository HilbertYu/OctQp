#ifndef RNDATA_H
#define RNDATA_H

#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
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

        typename std::vector<RnData<Dim> >::iterator begin(void)
        {
            return data_list.begin();
        }

        typename std::vector<RnData<Dim> >::iterator end(void)
        {
            return data_list.end();
        }

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

        static RnDataList<2> R2DataFileLoader(const std::string & file_name)
        {
            using namespace std;

            RnDataList<2> ret;

            ifstream ifs(file_name);
            if (ifs.bad())
            {
                fprintf(stderr, "error\n");
                exit(-1);
            }

            string line;
            while (getline(ifs, line))
            {
                if (ifs.eof())
                    break;

                if (line.size() == 0)
                    break;

                int coord[2] = {-1, -1};
                int r = sscanf(line.c_str(), "%d,%d\n", coord, coord+1);
                assert(r == 2);

                RnData<2> v;
                v(0) = coord[0];
                v(1) = coord[1];
                ret.push_back(v);
            }

            return ret;

        }
    };

};

#endif /* end of include guard: RNDATA_H */

