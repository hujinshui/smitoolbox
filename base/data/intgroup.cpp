/********************************************************************
 *
 *  intgroup.cpp
 *
 *  The C++ mex implementation of intgroup
 *
 *  Created by Dahua Lin, on June 6, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <vector>
#include <utility>

using namespace bcs;
using namespace bcs::matlab;


template<typename T>
inline void get_vs(const_marray mRgn, int& v0, int& v1)
{
    size_t ne = mRgn.nelems();
    const T *v = mRgn.data<T>();
        
    if (ne == 1)
    {
        v0 = 1;
        v1 = (int)(*v);
    }
    else
    {
        v0 = (int)(v[0]);
        v1 = (int)(v[1]);
    }
}



inline void get_range(const_marray mRgn, int& v0, int& v1)
{
    size_t ne = mRgn.nelems();
    
    if (!( (ne == 1 || ne == 2) && !mRgn.is_sparse() ) )
    {
        throw mexception("intgroup:invalidarg", 
            "The range [v0 v1] should be a (non-sparse) scalar or pair.");
    }
    
    if (mRgn.is_double())
    {
        get_vs<double>(mRgn, v0, v1);
    }
    else if (mRgn.is_single())
    {
        get_vs<float>(mRgn, v0, v1);
    }
    else if (mRgn.is_int32())
    {
        get_vs<int32_t>(mRgn, v0, v1);
    }
    else
    {
        throw mexception("intgroup:invalidarg", 
            "The range [v0 v1] should be of class double, single, or int32");
    }
}


template<typename T>
void group(int v0, int v1, const T *v, size_t n, std::vector<double>* gs)
{
    for (size_t i = 0; i < n; ++i)
    {
        int cv = (int)(v[i]);
        if (cv >= v0 && cv <= v1)
        {
            gs[cv - v0].push_back(i + 1);
        }
    }
}

inline void do_group(int v0, int v1, const_marray mVals, std::vector<double>* gs)
{        
    
    size_t n = mVals.nelems();
    
    if (mVals.is_double())
    {
        group(v0, v1, mVals.data<double>(), n, gs);
    }
    else if (mVals.is_single())
    {
        group(v0, v1, mVals.data<float>(), n, gs);
    }
    else if (mVals.is_int32())
    {
        group(v0, v1, mVals.data<int32_t>(), n, gs);
    }
    else
    {
        mexErrMsgIdAndTxt("intgroup:invalidarg", 
                "The class of values should be either double, single, int32");
    }    
}


marray groups_to_cells(size_t ng, const std::vector<double>* gs)
{    
    marray mC = create_mcell_array(1, ng);
    
    for (int i = 0; i < (int)ng; ++i)
    {
        mC.set_cell(i, to_matlab_row(gs[i]));
    }
    
    return mC;    
}


/**
 * Main entry:
 *
 *   Input:
 *     [0]: rgn:    the integer range in form of [v0, v1] (pair) 
 *     [1]: V:      the values array (double|single|int32)
 *   Output:
 *     [0]: G:      the cell array of indices for each value
 *
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
    {
        throw mexception("intgroup:invalidarg", 
            "The number of inputs to intgroup should be 2.");
    }

    const_marray mRgn(prhs[0]);
    const_marray mVals(prhs[1]);

    int v0, v1;
    get_range(mRgn, v0, v1);

    size_t ng = v1 - v0 + 1;
    std::vector<double> *vecs = new std::vector<double>[ng];

    do_group(v0, v1, mVals, vecs);

    plhs[0] = groups_to_cells(ng, vecs).mx_ptr();

    delete[] vecs;
}

BCSMEX_MAINDEF


