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
#include <bcslib/veccomp/vecstat.h>

using namespace bcs;
using namespace bcs::matlab;


template<typename T>
marray do_vec_median(const_marray mX)
{
    size_t n = mX.nelems();
    
    array1d<T> temp(n);
    temp << view1d<T>(mX);
    T v = median_inplace(n, temp.pbase());
    return create_mscalar<T>(v);
}


template<typename T>
marray do_mat_median(const_marray mX, int adim)
{
    size_t m = mX.nrows();
    size_t n = mX.ncolumns();
    
    const_aview2d<T, column_major_t> X = view2d<T>(mX);
    
    if (adim == 1)  // columns
    {
        marray mR = create_marray<T>(1, n);
        aview1d<T> r = view1d<T>(mR);
        
        array1d<T> temp(m);
        for (index_t i = 0; i < (index_t)n; ++i)
        {
            temp << X.column(i);
            r(i) = median_inplace(m, temp.pbase());            
        }
        
        return mR;        
    }
    else if (adim == 2)  // rows
    {
        marray mR = create_marray<T>(m, 1);
        aview1d<T> r = view1d<T>(mR);
        
        array1d<T> temp(n);
        for (index_t i = 0; i < (index_t)m; ++i)
        {
            temp << X.row(i);
            r(i) = median_inplace(n, temp.pbase());
        }
        
        return mR;
    }
}



template<typename T>
inline marray do_median(const_marray mX, int adim)
{
    size_t m = mX.nrows();
    size_t n = mX.ncolumns();
    
    if (adim == 0)
    {
        if (m == 1 || n == 1)
        {
            return do_vec_median<T>(mX);
        }
        else
        {
            return do_mat_median<T>(mX, 1);
        }            
    }
    
    
    if (adim == 1)
    {
        if (m == 1)
        {
            return duplicate(mX);            
        }
        else if (n == 1)
        {
            return do_vec_median<T>(mX);
        }
        else     
        {
            return do_mat_median<T>(mX, 1);
        }
    }
    else if (adim == 2)
    {
        if (n == 1)
        {
            return duplicate(mX);            
        }
        else if (m == 1)
        {
            return do_vec_median<T>(mX);
        }
        else
        {
            return do_mat_median<T>(mX, 2);
        }
    }
}



/**
 * main entry:
 *  [0] X:      the input matrix
 *  [1] dim:    the dimension along which median is solved
 *
 * output entry:
 *  [0] r:  the median values
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const_marray mX(prhs[0]);
    const_marray mDim(prhs[1]);
    
    int adim = (int)mDim.get_scalar<double>();
    
    if (mX.is_double())
    {
        plhs[0] = do_median<double>(mX, adim).mx_ptr();
    }
    else if (mX.is_single())
    {
        plhs[0] = do_median<float>(mX, adim).mx_ptr();
    }
}




