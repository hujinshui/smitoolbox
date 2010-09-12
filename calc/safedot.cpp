/********************************************************************
 *
 *  safedot.cpp
 *
 *  The C++ mex implementation of safe dot
 *
 ********************************************************************/

#include <mex.h>


template<typename T>
inline T safedot(int n, const T *x, const T* y)
{
    T s = T(0);
    for (int i = 0; i < n; ++i)
    {
        T cx = x[i];
        T cy = y[i];
        if (cx != 0 && cy != 0) s += cx * cy;
    }
    
    return s;
}

template<typename T>
inline T safedot(int n, const T *x, const T* y, int inc)
{
    T s = T(0);
    for (int i = 0; i < n; ++i)
    {
        T cx = *x;
        T cy = *y;
        if (cx != 0 && cy != 0) s += cx * cy;
        
        x += inc;
        y += inc;
    }
    
    return s;
}

template<typename T>
void safedot_cols(int m, int n, const T *A, const T *B, T *v)
{
    const T *a = A;
    const T *b = B;
    
    for (int i = 0; i < n; ++i)
    {
        v[i] = safedot(m, a, b);
        a += m;
        b += m;
    }
}


template<typename T>
void safedot_rows(int m, int n, const T *A, const T *B, T *v)
{
    const T *a = A;
    const T *b = B;
    
    for (int i = 0; i < m; ++i)
    {
        v[i] = safedot(n, a, b, m);
        ++a;
        ++b;
    }        
}


template<typename T>
mxArray* do_safedot(const mxArray *mxA, const mxArray *mxB, 
        mxClassID cid, int m, int n, int adim)
{
    const T *A = (const T*)mxGetData(mxA);
    const T *B = (const T*)mxGetData(mxB);
    
    mxArray *mxV = 0;
    
    if (adim == 0)
    {
        int len = n == 1 ? m : n;
        
        T v = safedot(len, A, B);
        mxV = mxCreateNumericMatrix(1, 1, cid, mxREAL);
        *((T*)mxGetData(mxV)) = v;
    }
    else if (adim == 1)
    {
        mxV = mxCreateNumericMatrix(1, n, cid, mxREAL);
        T *v = (T*)mxGetData(mxV);
        
        safedot_cols(m, n, A, B, v);
    }
    else if (adim == 2)
    {
        mxV = mxCreateNumericMatrix(m, 1, cid, mxREAL);
        T *v = (T*)mxGetData(mxV);
        
        safedot_rows(m, n, A, B, v);
    }
    
    return mxV;
}



inline bool is_real_matrix(const mxArray *mxA)
{
    return mxGetNumberOfDimensions(mxA) == 2 && !mxIsSparse(mxA) &&
            !mxIsComplex(mxA);
}

inline bool is_real_scalar(const mxArray *mxV)
{
    return mxGetNumberOfElements(mxV) == 1 && !mxIsComplex(mxV) && mxIsDouble(mxV);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // verify input
    
    if (nrhs < 2 || nrhs > 3)
    {
        mexErrMsgIdAndTxt("safedot:invalidarg", 
                "The number of input arguments to safedot is invalid.");
    }
    
    const mxArray *mxA = prhs[0];
    const mxArray *mxB = prhs[1];
    
    if (!is_real_matrix(mxA))
        mexErrMsgIdAndTxt("safedot:invalidarg", 
                "A should be a non-sparse real matrix.");
    
    if (!is_real_matrix(mxB))
        mexErrMsgIdAndTxt("safedot:invalidarg",
                "B should be a non-sparse real matrix.");
    
    int m = (int)mxGetM(mxA);
    int n = (int)mxGetN(mxA);
    
    if (mxGetM(mxB) != m || mxGetN(mxB) != n)
    {
        mexErrMsgIdAndTxt("safedot:invalidarg", 
                "The sizes of A and B are inconsistent.");
    }
    
    // decide along dimension
    
    int along_dim = 0;
    if (nrhs < 3)
    {
        if (m > 1) along_dim = 1;
    }    
    else
    {
        const mxArray *mxDim = prhs[2];
        if (!is_real_scalar(mxDim))
        {
            mexErrMsgIdAndTxt("safedot:invalidarg", "dim should be an integer scalar.");
        }
        
        int dim = (int)mxGetScalar(mxDim);
        if (dim < 1 || dim > 2)
        {
            mexErrMsgIdAndTxt("safedot:invalidarg", "dim should be either 1 or 2.");
        }
        
        if (dim == 1)
        {
            along_dim = n == 1 ? 0 : dim;
        }
        else
        {
            along_dim = m == 1 ? 0 : dim;
        }
    }
    
    // do computation
    
    if (mxIsDouble(mxA) && mxIsDouble(mxB))
    {
        plhs[0] = do_safedot<double>(mxA, mxB, mxDOUBLE_CLASS, m, n, along_dim);
    }
    else if (mxIsSingle(mxA) && mxIsSingle(mxB))
    {
        plhs[0] = do_safedot<double>(mxA, mxB, mxSINGLE_CLASS, m, n, along_dim);
    }
    else
    {
        mexErrMsgIdAndTxt("safedot:invalidarg", 
                "A and B should be both double or both single.");
    }
}




