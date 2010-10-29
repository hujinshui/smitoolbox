/********************************************************************
 *
 *  safedot.cpp
 *
 *  The C++ mex implementation of safe dot
 *
 ********************************************************************/


#include "../clib/marray.h"

using namespace smi;

template<typename T>
inline T safedot(const VectorCView<T>& veca, const VectorCView<T>& vecb)
{
    int n = veca.nelems();
    
    const T *a = veca.data();
    const T *b = vecb.data();
    
    int iv_a = veca.interval();
    int iv_b = vecb.interval();
    
    T s = T(0);
    for (int i = 0; i < n; ++i)
    {
        T av = *a;
        T bv = *b;
        
        if (av != 0 && bv != 0)
        {
            s += av * bv;
        }
        
        a += iv_a;
        b += iv_b;
    }
    
    return s;
}


template<typename T>
inline mxArray* create_vector(int n, int dim)
{
    return dim == 1 ? create_matlab_matrix<T>(1, n) : create_matlab_matrix<T>(n, 1);
}


template<typename T>
mxArray* multi_safedot(const MultiVectorCView<T>& A, const MultiVectorCView<T>& B, int dim)
{
    if (A.veclen() != B.veclen())
        mexErrMsgIdAndTxt("safedot:invalidarg", "The vector dimensions are inconsistent.");
    
    int na = A.nvecs();
    int nb = B.nvecs();
    
    mxArray *mxR = 0;
    
    if (na == 1)
    {
        if (nb == 1)
        {
            T r = safedot(A[0], B[0]);
            mxR = create_matlab_scalar<T>(r);
        }
        else
        {
            mxR = create_vector<T>(nb, dim);
            VectorView<T> R = MArray(mxR).to_vector<T>();
            
            VectorCView<T> a = A[0];
            
            for (int i = 0; i < nb; ++i)            
                R(i) = safedot(a, B[i]);            
        }
    }
    else
    {
        if (nb == 1)
        {
            mxR = create_vector<T>(na, dim);
            VectorView<T> R = MArray(mxR).to_vector<T>();
            
            VectorCView<T> b = B[0];
            
            for (int i = 0; i < na; ++i)
                R(i) = safedot(A[i], b);            
        }
        else if (na == nb)
        {
            mxR = create_vector<T>(na, dim);
            VectorView<T> R = MArray(mxR).to_vector<T>();
            
            for (int i = 0; i < na; ++i)
                R(i) = safedot(A[i], B[i]);
        }
        else
        {
            mexErrMsgIdAndTxt("safedot:invalidarg", 
                    "The sizes of A and B are inconsistent.");
        }        
    }
    
    return mxR;
}

template<typename T>
mxArray* run(const MArray& mA, const MArray& mB, int dim)
{
    MatrixCView<T> A = mA.to_matrix<T>();        
    MatrixCView<T> B = mB.to_matrix<T>();
    
    if (dim == 1)
    {
        return multi_safedot<T>(A.columns(), B.columns(), dim);
    }
    else
    {
        return multi_safedot<T>(A.rows(), B.rows(), dim);
    }
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    if (nrhs < 2 || nrhs > 3)
        mexErrMsgIdAndTxt("safedot:invalidarg", 
                "The number of inputs should be either 2 or 3.");
    
    MArray mA(prhs[0]);
    MArray mB(prhs[1]);
    
    if (!mA.is_real_nonsparse_matrix())
        mexErrMsgIdAndTxt("safedot:invalidarg", "A should be a real non-sparse matrix.");
    
    if (!mB.is_real_nonsparse_matrix())
        mexErrMsgIdAndTxt("safedot:invalidarg", "B should be a real non-sparse matrix.");
    
    int dim = 0;
    if (nrhs >= 3)
    {
        MArray mDim(prhs[2]);
        if (!mDim.is_real_double_scalar())
            mexErrMsgIdAndTxt("safedot:invalidarg", "dim should be a double scalar.");
        
        dim = (int)mDim.get_double_scalar();
        if (!(dim == 1 || dim == 2))
            mexErrMsgIdAndTxt("safedot:invalidarg", "dim should be either 1 or 2.");                
    }
    
    // determine dim
    
    if (dim == 0)
    {
        dim = (mA.nrows() == 1 || mB.nrows() == 1) ? 2 : 1;           
    }
    
    // main
    
    if (mA.is_double() && mB.is_double())
    {
        plhs[0] = run<double>(mA, mB, dim);                
    }
    else if (mA.is_single() && mB.is_single())
    {
        plhs[0] = run<float>(mA, mB, dim);
    }
    else
    {
        mexErrMsgIdAndTxt("safedot:invalidarg", 
                "A and B should be both double or both single.");
    }    
    
}




