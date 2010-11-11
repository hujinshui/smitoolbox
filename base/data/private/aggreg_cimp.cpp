/********************************************************************
 *
 *  aggreg_cimp.cpp
 *
 *  The C++ mex implementation of aggreg.m
 *
 *  Created by Dahua Lin, on Nov 10, 2010
 *
 ********************************************************************/

#include <mex.h>
#include <limits>
#include <algorithm>

// recursive functor for aggregation

template<typename T>
struct rsum
{
    inline void operator()(T& s, T v) const
    {
        s += v;
    }
};

template<typename T>
struct rmin
{   
    inline void operator()(T& s, T v) const
    {
        if (v < s) s = v;
    }
};

template<typename T>
struct rmax
{
    inline void operator()(T& s, T v) const
    {   
        if (v > s) s = v;
    }
};


// aggregation functions


template<typename T, typename RFun>
void aggregate_rows(int m, int n, int K, const T *X, const int *I, T *R)
{
    RFun rfun;
    
    for (int i = 0; i < m; ++i)
    {
        int k = I[i];
        if (k >= 0 && k < K)
        {
            const T *x = X + i;
            T *r = R + k;
            
            for (int j = 0; j < n; ++j)
            {
                rfun(*r, *x);
                
                x += m;
                r += K;
            }
        }
    }
}


template<typename T, typename RFun>
void aggregate_cols(int m, int n, int K, const T *X, const int *I, T *R)
{
    RFun rfun;
    
    for (int i = 0; i < n; ++i)
    {
        int k = I[i];
        if (k >= 0 && k < K)
        {
            const T *x = X + i * m;
            T *r = R + k * m;
            
            for (int j = 0; j < m; ++j)
            {
                rfun(r[j], x[j]);
            }
        }
    }
}


template<typename T>
void aggregate(int m, int n, int K, const T *X, const int *I, T *R,
        bool byrows, int code)
{
    if (code == 1)
    {
        if (byrows)
            aggregate_rows<T, rsum<T> >(m, n, K, X, I, R);
        else
            aggregate_cols<T, rsum<T> >(m, n, K, X, I, R);
    }
    else if (code == 2)
    {
        if (byrows)
            aggregate_rows<T, rmin<T> >(m, n, K, X, I, R);
        else
            aggregate_cols<T, rmin<T> >(m, n, K, X, I, R);
    }
    else if (code == 3)
    {
        if (byrows)
            aggregate_rows<T, rmax<T> >(m, n, K, X, I, R);
        else
            aggregate_cols<T, rmax<T> >(m, n, K, X, I, R);
    }
}


template<typename T>
inline T get_initv(int code)
{
    if (code == 1) // sum
    {
        return T(0);
    }
    else if (code == 2) // min
    {
        return std::numeric_limits<T>::infinity();
    }
    else // max
    {
        return - std::numeric_limits<T>::infinity();
    }
}

template<typename T> mxArray* init_matrix(int m, int n, int code);  

template<>
inline mxArray* init_matrix<double>(int m, int n, int code)
{
    double v0 = get_initv<double>(code);
    mxArray *mx = mxCreateDoubleMatrix(m, n, mxREAL);
    std::fill_n(mxGetPr(mx), m*n, v0);
    return mx;
}

template<>
inline mxArray* init_matrix<float>(int m, int n, int code)
{
    float v0 = get_initv<float>(code);
    mxArray *mx = mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
    std::fill_n((float*)mxGetData(mx), m*n, v0);
    return mx;
}




/**
 * main entry:
 *
 * Inputs:
 *  [0]:  X:    data
 *  [1]:  K:    #classes
 *  [2]:  I:    indices
 *  [3]:  code: function code (1-sum, 2-min, 3-max)
 *
 * Outputs:
 *  [0]:  R:    the results
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxX = prhs[0];
    const mxArray *mxK = prhs[1];
    const mxArray *mxI = prhs[2];
    const mxArray *mxCode = prhs[3];
    
    int m = (int)mxGetM(mxX);
    int n = (int)mxGetN(mxX);
    
    bool byrows = mxGetM(mxI) != 1;    
    int K = (int)mxGetScalar(mxK);
    int code = (int)mxGetScalar(mxCode);
    
    // prepare output
    
    mxArray *mxR = 0;
    
    int mr = byrows ? K : m;
    int nr = byrows ? n : K;
    
    // main
    
    const int *I = (const int*)mxGetData(mxI);
    
    switch (mxGetClassID(mxX))
    {
        case mxDOUBLE_CLASS:
            mxR = init_matrix<double>(mr, nr, code);
            aggregate(m, n, K, (const double*)mxGetData(mxX), I, 
                    (double*)mxGetData(mxR), byrows, code);
            break;
            
        case mxSINGLE_CLASS:
            mxR = init_matrix<float>(mr, nr, code);
            aggregate(m, n, K, (const float*)mxGetData(mxX), I, 
                    (float*)mxGetData(mxR), byrows, code);
            break;   
            
        default:
            mexErrMsgIdAndTxt("aggreg:invalidarg", 
                    "The type of X must be either double or single.");
    }
    
    plhs[0] = mxR;
    
}


