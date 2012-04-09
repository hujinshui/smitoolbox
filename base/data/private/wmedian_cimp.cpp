/********************************************************************
 *
 *  wmedian_cimp.cpp
 *
 *  The implementation of the core part of wmedian
 *
 *  Created by Dahua Lin, on Sep 28, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

using namespace bcs;
using namespace bcs::matlab;


template<typename T>
inline T find_wmedian(
        const caview1d<T>& x,   // x : sorted values
        const caview1d<T>& f)   // f : cumsum of sorted weights
{
    int n = (int)x.nelems();
    T t = f[n-1] / 2;
    
    int i = 0;
    while (i < n && f[i] < t) ++i;
    
    if (f[i] > t)
    {
        return x[i];
    }
    else
    {
        return i < n-1 ? (x[i] + x[i+1]) / 2 : x[i];
    }
}


template<typename T>
inline marray do_wmedian(const_marray mX, const_marray mF)
{
    size_t m = mX.nrows();
    size_t n = mX.ncolumns();
    
    if (m > 1)
    {
        if (n == 1)
        {
            caview1d<T> x = view1d<T>(mX);
            caview1d<T> f = view1d<T>(mF);
            
            T r = find_wmedian(x, f);
            
            return create_mscalar<T>(r);
        }
        else
        {
            caview2d<T> X = view2d<T>(mX);
            caview2d<T> F = view2d<T>(mF);
            
            marray mR = create_marray<T>(1, n);
            aview1d<T> R = view1d<T>(mR);
            
            for (index_t i = 0; i < (index_t)n; ++i)
            {
                R(i) = find_wmedian(X.column(i), F.column(i));
            }
            
            return mR;
        }
    }
    else
    {
        return duplicate(mX);
    }
}



/**
 * main entry
 * 
 * Input
 *  [0]:     X   (sorted values) (along dim 1)
 *  [1]:     F   (cumsum of sorted weights)
 * Output
 *  [0]:     R   (weighted median values)
 *
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const_marray mX(prhs[0]);
    const_marray mF(prhs[1]);
    
    marray mR;
    if (mX.is_double())
    {
        mR = do_wmedian<double>(mX, mF);
    }
    else
    {
        mR = do_wmedian<float>(mX, mF);
    }
    
    plhs[0] = mR.mx_ptr();
}


BCSMEX_MAINDEF





