/********************************************************************
 *
 *  fast_median_cimp.cpp
 *
 *  The C++ mex implementation for fast_median.m
 *
 *  Created by Dahua Lin, on Nov 11, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/array/array2d.h>
#include <algorithm>

using namespace bcs;
using namespace bcs::matlab;

template<typename T>
inline T my_median_inplace(index_t n, T *x)
{    
    if (n == 1)
    {
        return *x;
    }
    else if (n == 2)
    {
        T x0 = x[0];
        T x1 = x[1];
        return x0 + (x1 - x0) / 2;
    }
    else if (n % 2 == 0) // even
    {
        T *pm = x + (n/2);
        std::nth_element(x, pm, x+n);
        
        T v1 = *pm;
        T v0 = *(std::max_element(x, pm));
        
        return v0 + (v1 - v0) / 2;
    }
    else  // odd
    {
        T *pm = x + (n/2);
        std::nth_element(x, pm, x+n);
        return *pm;
    }
}


template<typename T>
inline marray do_median(const_marray mX, int d)
{
    if (d == 0)  // for vec
    {
        index_t n = (index_t)mX.nelems();
        array1d<T> temp(n, mX.data<T>());
        T v = my_median_inplace(n, temp.pbase());
        
        return create_mscalar<T>(v);
    }
    else // for mat
    {
        index_t m = (index_t)mX.nrows();
        index_t n = (index_t)mX.ncolumns();
        
        array2d<T> temp(m, n, mX.data<T>());
        
        marray mR = create_marray<T>(1, (size_t)n);
        aview1d<T> r = view1d<T>(mR);
        
        for (index_t i = 0; i < n; ++i)
        {
            r(i) = my_median_inplace(m, &(temp(0, i)));            
        }
        
        return mR;
    }
}



/**
 * main entry:
 *  [0] X:      the input matrix
 *  [1] d:      d == 0: solve median of a vector
 *              d == 1: solve median(s) of each column
 *
 * output entry:
 *  [0] r:  the median values
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const_marray mX(prhs[0]);
    const_marray mDim(prhs[1]);
    
    int d = (int)mDim.get_scalar<double>();
    
    if (mX.is_double())
    {
        plhs[0] = do_median<double>(mX, d).mx_ptr();
    }
    else if (mX.is_single())
    {
        plhs[0] = do_median<float>(mX, d).mx_ptr();
    }
}




