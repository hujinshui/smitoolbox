/**********************************************************
 *
 *  sum_xlogy.cpp
 *
 *  The C++ mex implementation of sum_xlogy
 *
 *  Created by Dahua Lin, on Mar 28, 2011
 *
 **********************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <cmath>

using namespace bcs;
using namespace bcs::matlab;


inline int get_dim(const_marray mDim)
{
    if ( !(mDim.is_scalar() && mDim.is_numeric()) )
    {
        throw mexception("sum_xlogy:invalidarg", 
            "dim should be a numeric scalar.");
    }
    
    int d = 0;
    switch (mDim.class_id())
    {
        case mxDOUBLE_CLASS:
            d = (int)mDim.get_scalar<double>();
            break;
        case mxSINGLE_CLASS:
            d = (int)mDim.get_scalar<float>();
            break;
        case mxINT32_CLASS:
            d = (int)mDim.get_scalar<int32_t>();
            break;
        case mxUINT32_CLASS:
            d = (int)mDim.get_scalar<uint32_t>();
            break;
        default:
            throw mexception("sum_xlogy:invalidarg", 
                "dim should be of type double, single, int32, or uint32.");            
    }
    
    if ( !(d == 1 || d == 2) )
    {
        throw mexception("sum_xlogy:invalidarg", 
            "dim should be either 1 or 2.");                
    }
    
    return d;
}


template<typename T, class TIndexerX, class TIndexerY>
T sum_xlogy(const const_aview1d<T, TIndexerX>& x, const const_aview1d<T, TIndexerY>& y)
{
    T s(0);
    index_t n = (index_t)x.nelems();
    
    for (index_t i = 0; i < n; ++i)
    {
        T xi = x(i);
        if (xi != 0)
        {
            s += xi * std::log(y(i));
        }
    }
        
    return s;
}


template<typename T>
inline marray do_sum_xlogy_columns(const_marray mX, const_marray mY)
{
    size_t x_n = mX.ncolumns();
    size_t y_n = mY.ncolumns();
    
    size_t n = max(x_n, y_n);
    
    if (n == 1)
    {
        const_aview1d<T> x = view1d<T>(mX);        
        const_aview1d<T> y = view1d<T>(mY);        
        T r = sum_xlogy(x, y);
        return create_mscalar<T>(r);
    }
    else
    {
        marray mR = create_marray<T>(1, n);
        aview1d<T> r = view1d<T>(mR);
        
        if (x_n == 1)
        {
            const_aview1d<T> x = view1d<T>(mX); 
            const_aview2d<T, column_major_t> Y = view2d<T>(mY);
            
            for (index_t i = 0; i < (index_t)n; ++i)
            {
                r(i) = sum_xlogy(x, Y.column(i));
            }            
        }
        else if (y_n == 1)
        {
            const_aview2d<T, column_major_t> X = view2d<T>(mX);
            const_aview1d<T> y = view1d<T>(mY); 
            
            for (index_t i = 0; i < (index_t)n; ++i)
            {
                r(i) = sum_xlogy(X.column(i), y);
            }
        }
        else
        {
            const_aview2d<T, column_major_t> X = view2d<T>(mX);
            const_aview2d<T, column_major_t> Y = view2d<T>(mY);
            
            for (index_t i = 0; i < (index_t)n; ++i)
            {
                r(i) = sum_xlogy(X.column(i), Y.column(i));
            }
        }
        
        return mR;
    }  
}


template<typename T>
inline marray do_sum_xlogy_rows(const_marray mX, const_marray mY)
{
    size_t x_m = mX.nrows();
    size_t y_m = mY.nrows();
    
    size_t m = max(x_m, y_m);
    
    if (m == 1)
    {
        const_aview1d<T> x = view1d<T>(mX);        
        const_aview1d<T> y = view1d<T>(mY);        
        T r = sum_xlogy(x, y);
        return create_mscalar<T>(r);
    }
    else
    {
        marray mR = create_marray<T>(m, 1);
        aview1d<T> r = view1d<T>(mR);
        
        if (x_m == 1)
        {
            const_aview1d<T> x = view1d<T>(mX); 
            const_aview2d<T, column_major_t> Y = view2d<T>(mY);
            
            for (index_t i = 0; i < (index_t)m; ++i)
            {
                r(i) = sum_xlogy(x, Y.row(i));
            }            
        }
        else if (y_m == 1)
        {
            const_aview2d<T, column_major_t> X = view2d<T>(mX);
            const_aview1d<T> y = view1d<T>(mY); 
            
            for (index_t i = 0; i < (index_t)m; ++i)
            {
                r(i) = sum_xlogy(X.row(i), y);
            }
        }
        else
        {
            const_aview2d<T, column_major_t> X = view2d<T>(mX);
            const_aview2d<T, column_major_t> Y = view2d<T>(mY);
            
            for (index_t i = 0; i < (index_t)m; ++i)
            {
                r(i) = sum_xlogy(X.row(i), Y.row(i));
            }
        }
        
        return mR;
    }         
}



/**
 * main entry:
 *
 * Input
 *   [0] X:    the matrix of x values
 *   [1] Y:    the matrix of y values
 *   [2] dim:  the dimension along which the sum is taken (optional)
 *
 * Output
 *   [0] R:    the result matrix
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    if (nrhs < 2 || nrhs > 3)
    {
        throw mexception("sum_xlogy:invalidarg", 
            "The number of inputs to sum_xlogy should be 2 or 3.");
    }
    
    const_marray mX(prhs[0]);
    const_marray mY(prhs[1]);
    
    if ( !( mX.is_matrix() && mX.is_float() && !mX.is_sparse() && !mX.is_complex() ))
    {
        throw mexception("sum_xlogy:invalidarg", 
            "X should be a float-typed non-sparse non-complex matrix.");
    }
    if ( !( mY.is_matrix() && mY.is_float() && !mY.is_sparse() && !mY.is_complex() ))
    {
        throw mexception("sum_xlogy:invalidarg", 
            "Y should be a float-typed non-sparse non-complex matrix.");
    }        
    
    if ( mX.class_id() != mY.class_id() )
    {
        throw mexception("sum_xlogy:invalidarg", 
            "X and Y should have the same value types.");
    }        
    
    // determine adim
    
    size_t x_m = mX.nrows();
    size_t x_n = mX.ncolumns();
    size_t y_m = mY.nrows();
    size_t y_n = mY.ncolumns();
    
    int adim = 0;    
    if (nrhs >= 3)
    {
        const_marray mDim(prhs[2]);
        adim = get_dim(mDim);      
        
        if (adim == 1)
        {
            if (!( (x_m == y_m) && (x_n == y_n || x_n == 1 || y_n == 1) ))
            {
                adim = 0; // tag as invalid
            }
        }
        else // adim == 2
        {
            if (!( (x_n == y_n) && (x_m == y_m || x_m == 1 || y_m == 1) ))
            {
                adim = 0;
            }
        }
    }
    else
    {
        if (x_m == y_m && x_n == y_n)
        {
            adim = (x_m == 1 ? 2 : 1);
        }
        else
        {
            if (x_m == y_m && (x_n == 1 || y_n == 1))
            {
                adim = 1;
            }
            else if (x_n == y_n && (x_m == 1 || y_m == 1))
            {
                adim = 2;
            }
        }
    }
    
    if (adim == 0)
    {
        throw mexception("sum_xlogy:invalidarg", 
            "The sizes of X and Y are inconsistent.");
    }
    
    
    // main delegation
    
    marray mR;
    
    if (adim == 1)
    {
        if (mX.is_double())
        {
            mR = do_sum_xlogy_columns<double>(mX, mY);
        }
        else
        {
            mR = do_sum_xlogy_columns<float>(mX, mY);
        }        
    }
    else
    {
        if (mX.is_double())
        {
            mR = do_sum_xlogy_rows<double>(mX, mY);
        }
        else
        {
            mR = do_sum_xlogy_rows<float>(mX, mY);
        }         
    }
    
    // output
    
    plhs[0] = mR.mx_ptr();    
}

BCSMEX_MAINDEF

        
        
