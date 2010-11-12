/********************************************************************
 *
 *  safedot.cpp
 *
 *  The C++ mex implementation of safe dot
 *
 ********************************************************************/

#include <mex.h>

template<typename T>
inline T smul(T x, T y)
{
    return (x == 0 || y == 0) ? 0 : x * y;
}



template<typename T>
inline T vsafedot(int n, const T *a, const T *b, int intv)
{
    T s(0);
    
    if (intv == 1)
    {
        for (int i = 0; i < n; ++i)
        {
            s += smul(a[i], b[i]);
        }
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            s += smul((*a), (*b));
            a += intv;
            b += intv;
        }
    }
    
    return s;
}


template<typename T>
inline void safedot_cols(int m, int n, const T *a, const T *b, T *r)
{
    for (int i = 0; i < n; ++i)
    {
        r[i] = vsafedot(m, a+i*m, b+i*m, 1);
    }
}


template<typename T>
inline void safedot_rows(int m, int n, const T *a, const T *b, T *r)
{
    for (int i = 0; i < m; ++i)
    {
        r[i] = vsafedot(n, a+i, b+i, m);
    }
}


inline int decide_dim(int m, int n, int dim)
{
    if (dim == 0)
    {
        return (m == 1 || n == 1) ? 0 : 1;
    }
    else
    {
        return ((dim == 1 && n == 1) || (dim == 2 && m == 1)) ? 0 : dim;
    }
}




template<typename T>
inline mxArray* do_safedot(const mxArray *mxA, const mxArray *mxB, 
        mxClassID cid, int m, int n, int dim)
{    
    int d = decide_dim(m, n, dim);
    
    mxArray *mxR = 0;    
    
    const T *a = (const T*)mxGetData(mxA);
    const T *b = (const T*)mxGetData(mxB);
    
    
    if (d == 0)
    {
        mxR = mxCreateNumericMatrix(1, 1, cid, mxREAL);
        *((T*)mxGetData(mxR)) = vsafedot(m*n, a, b, 1);
    }
    else if (d == 1)
    {
        mxR = mxCreateNumericMatrix(1, n, cid, mxREAL);
        T *r = (T*)mxGetData(mxR);
        safedot_cols(m, n, a, b, r);
    }
    else
    {
        mxR = mxCreateNumericMatrix(m, 1, cid, mxREAL);
        T *r = (T*)mxGetData(mxR);
        safedot_rows(m, n, a, b, r);
    }
    
    return mxR;
    
}


inline bool is_nonsparse_realmat(const mxArray *mxA)
{
    return mxIsNumeric(mxA) && 
            mxGetNumberOfDimensions(mxA) == 2 &&
            !mxIsSparse(mxA) &&
            !mxIsComplex(mxA);            
}





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take and verify input
    
    if (nrhs < 2 || nrhs > 3)
        mexErrMsgIdAndTxt("safedot:invalidarg", 
                "The number of inputs should be either 2 or 3.");
    
    const mxArray *mxA = prhs[0];
    const mxArray *mxB = prhs[1];
    
    if (!is_nonsparse_realmat(mxA))
        mexErrMsgIdAndTxt("safedot:invalidarg", "A should be a real non-sparse matrix.");
    
    if (!is_nonsparse_realmat(mxB))
        mexErrMsgIdAndTxt("safedot:invalidarg", "B should be a real non-sparse matrix.");
    
    int m1 = mxGetM(mxA);
    int n1 = mxGetN(mxA);
    int m2 = mxGetM(mxB);
    int n2 = mxGetN(mxB);       
    
    if (!(m1 == m2 && n1 == n2))
        mexErrMsgIdAndTxt("safedot:invalidarg", "The sizes of A and B are inconsistent.");
    
    int m = m1;
    int n = n1;
    
    int dim = 0;
    if (nrhs >= 3)
    {
        const mxArray *mxDim = prhs[2];
        if ( !(mxIsDouble(mxDim) && mxGetNumberOfElements(mxDim) == 1) )
            mexErrMsgIdAndTxt("safedot:invalidarg", "dim should be a double scalar.");
        
        dim = (int)mxGetScalar(mxDim);
        if (!(dim == 1 || dim == 2))
            mexErrMsgIdAndTxt("safedot:invalidarg", "dim should be either 1 or 2.");                
    }
    
    // main
    
    mxClassID cid = mxGetClassID(mxA);
    mxClassID cid2 = mxGetClassID(mxB);
    
    if (cid != cid2)
        mexErrMsgIdAndTxt("safedot:invalidarg", "A and B should have the same value type.");
    
    mxArray *mxR = 0;
    
    switch (cid)
    {
        case mxDOUBLE_CLASS:
            mxR = do_safedot<double>(mxA, mxB, cid, m, n, dim);
            break;
            
        case mxSINGLE_CLASS:
            mxR = do_safedot<float>(mxA, mxB, cid, m, n, dim);
            break;
            
        default:
            mexErrMsgIdAndTxt("safedot:invalidarg", 
                    "A and B should be either double or single.");
    }
    
    plhs[0] = mxR;
    
}




