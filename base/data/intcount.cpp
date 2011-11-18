/********************************************************************
 *
 *  intcount.cpp
 *
 *  The C++ mex implementation of integer counting
 *
 *  Created by Dahua Lin, on May 26, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

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
        throw mexception("intcount:invalidarg", 
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
        throw mexception("intcount:invalidarg", 
            "The range [v0 v1] should be of class double, single, or int32");
    }
}


template<typename T>
void count(int v0, int v1, const T *v, size_t n, double *c)
{
    for (size_t i = 0; i < n; ++i)
    {
        int cv = (int)(v[i]);
        if (cv >= v0 && cv <= v1)
        {
            ++ c[cv - v0];
        }
    }
}


inline marray do_count(int v0, int v1, const_marray mVals)
{
    int m = v1 - v0 + 1;
    marray mCount = create_marray<double>(1, m);
    double *c = mCount.data<double>();
    
    size_t n = mVals.nelems();
    
    if (mVals.is_double())
    {
        count(v0, v1, mVals.data<double>(), n, c);
    }
    else if (mVals.is_single())
    {
        count(v0, v1, mVals.data<float>(), n, c);
    }
    else if (mVals.is_int32())
    {
        count(v0, v1, mVals.data<int32_t>(), n, c);
    }
    else
    {
        throw mexception("intcount:invalidarg", 
            "The class of values should be either double, single, int32");
    }
    
    return mCount;
    
}


/***********
 *
 *  Main entry
 *   
 *  Input:
 *    [0]:  the value range in form of [v0, v1] [pair (double|single|int32)]
 *    [1]:  the array of values to count [array (double|single|int32)]
 *
 *  Output:
 *    [0]:  the vector of counts (of size 1 x (v1 - v0 + 1)) [double array]
 *
 *  No 
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
    {
        throw mexception("intcount:invalidarg", 
            "The number of inputs to intcount should be 2.");
    }

    const_marray mRgn(prhs[0]);
    const_marray mVals(prhs[1]);

    int v0, v1;
    get_range(mRgn, v0, v1);

    plhs[0] = do_count(v0, v1, mVals).mx_ptr();

}


BCSMEX_MAINDEF


