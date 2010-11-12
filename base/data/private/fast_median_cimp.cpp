/********************************************************************
 *
 *  fast_median_cimp.cpp
 *
 *  The C++ mex implementation for fast_median.m
 *
 *  Created by Dahua Lin, on Nov 11, 2010
 *
 ********************************************************************/

#include <mex.h>
#include <string.h>
#include <algorithm>


template<typename T>
T two_mean(T x0, T x1)
{
    return x0 + (x1 - x0) / 2;
}


template<typename T>
inline T find_median_odd(T *x, int n)        
{    
    if (n == 1)
    {
        return *x;
    }
    else
    {
        T *pm = x + (n/2);
        std::nth_element(x, pm, x+n);
        return *pm;
    }
}

template<typename T>
inline T find_median_even(T *x, int n)
{
    if (n == 2)
    {
        return two_mean(x[0], x[1]);
    }
    else
    {
        T *pm = x + (n/2);
        std::nth_element(x, pm, x+n);
    
        T v1 = *pm;
        T v0 = *(std::max_element(x, pm));
    
        return two_mean(v0, v1);
    }
}


template<typename T>
inline T vec_median(const T *x, int n)
{        
    if (n == 1)
    {
        return *x;
    }
    else if (n == 2)
    {
        return two_mean(x[0], x[1]);
    }
    else
    {
        T *y = new T[n];
        ::memcpy(y, x, n * sizeof(T));
    
        T v = n % 2 == 0 ? 
            find_median_even(y, n) : 
            find_median_odd(y, n);
            
        delete[] y;            
        return v;
    }
}


template<typename T>
void mat_median(const T *x, int m, int n, T *r)
{
    T *y = new T[m * n];
    ::memcpy(y, x, m * n * sizeof(T));
    
    if (m % 2 == 0)
    {
        for (int i = 0; i < n; ++i)
        {
            r[i] = find_median_even(y+i*m, m);
        }
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            r[i] = find_median_odd(y+i*m, m);
        }
    }
    
    delete[] y;
}


template<typename T> mxArray* make_mat(int m, int n);

template<> mxArray* make_mat<double>(int m, int n)
{
    return mxCreateNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);
}

template<> mxArray* make_mat<float>(int m, int n)
{
    return mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
}



template<typename T>
mxArray* do_median(const mxArray *mxX)
{
    int m = (int)mxGetM(mxX);
    int n = (int)mxGetN(mxX);
    
    mxArray *mxR = 0;
    
    if (m == 1)
    {
        mxR = make_mat<T>(1, 1);
        *((T*)mxGetData(mxR)) = vec_median((const T*)mxGetData(mxX), n);
    }
    else if (n == 1)
    {
        mxR = make_mat<T>(1, 1);
        *((T*)mxGetData(mxR)) = vec_median((const T*)mxGetData(mxX), m);
    }
    else
    {
        mxR = make_mat<T>(1, n);
        mat_median((const T*)mxGetData(mxX), m, n, (T*)mxGetData(mxR));
    }
    
    return mxR;
}



/**
 * main entry:
 *  [0] X:  the input matrix
 *
 * output entry:
 *  [0] r:  the median values
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *mxX = prhs[0];
    
    if (mxIsDouble(mxX))
    {
        plhs[0] = do_median<double>(mxX);
    }
    else if (mxIsSingle(mxX))
    {
        plhs[0] = do_median<float>(mxX);
    }
}




