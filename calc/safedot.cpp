/********************************************************************
 *
 *  safedot.cpp
 *
 *  The C++ mex implementation of safe dot
 *
 ********************************************************************/

#include <mex.h>


struct VecSpec
{
    int num_vecs;   // the number of vectors
    int vlen;   // the vector length
    int intv;   // element interval
    int step;   // forward offset for next vector
    
    void set_rows(int m, int n, int num)
    {
        num_vecs = num;
        vlen = n;
        intv = m;
        step = m == 1 ? 0 : 1;
    }
    
    void set_cols(int m, int n, int num)
    {
        num_vecs = num;
        vlen = m;
        intv = 1;
        step = n == 1 ? 0 : m;
    }
};



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
    
    int dim = 0;
    if (nrhs >= 3)
    {
        const mxArray *mxDim = prhs[2];
        if (!is_real_scalar(mxDim))
        {
            mexErrMsgIdAndTxt("safedot:invalidarg", "dim should be an integer scalar.");
        }
        dim = (int)mxGetScalar(mxDim);
        
        if (dim < 1 || dim > 2)
        {
            mexErrMsgIdAndTxt("safedot:invalidarg", "dim should be either 1 or 2.");
        }
    }
    
    
    // determine dimensions
    
    int m1 = (int)mxGetM(mxA);
    int n1 = (int)mxGetN(mxA);
    
    int m2 = (int)mxGetM(mxB);
    int n2 = (int)mxGetN(mxB);
    
    VecSpec S1;
    VecSpec S2;
    
    if (dim == 0)
    {
        if ((m1 == m2
    }    
    
    if (dim == 1)
    {
        if (m1 != m2)
            mexErrMsgIdAndTxt("safedot:invalidarg", 
                    "A and B have different lengths along the specified dimension.");    
        
        if (!(n1 == n2 || n1 == 1 || n2 == 1))
            mexErrMsgIdAndTxt("safedot:invalidarg", 
                    "The number of vectors in A and B are inconsistent.");
        
        int num = n1 > n2 ? n1 : n2;
        
        S1.set_cols(m1, n1, num);
        S2.set_cols(m2, n2, num);
    }
    else // dim == 2
    {
        if (n1 != n2)
            mexErrMsgIdAndTxt("safedot:invalidarg",
                    "A and B have different lengths along the specified dimension."); 
        
        if (!(m1 == m2 || m1 == 1 || m2 == 1))
            mexErrMsgIdAndTxt("safedot:invalidarg", 
                    "The number of vectors in A and B are inconsistent.");
        
        int num = m1 > m2 ? m1 : m2;
        
        S1.set_rows(m2, n2, num);
        S2.set_rows(m2, n2, num);
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




