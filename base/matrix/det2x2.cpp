/********************************************************************
 *
 *  det2x2.cpp
 *
 *  The C++ mex implementation of determinant of 2x2 matrices
 *
 *  Created by Dahua Lin, on Jun 8, 2010
 *
 ********************************************************************/


#include "smallmat.h"


template<typename T>
inline mxArray* do_det2x2(const mxArray *mxA)
{    
    int n = 0;
    
    if (test_size(mxA, 2, 2, n))
    {        
        const Mat2x2<T>* A = (const Mat2x2<T>*)mxGetData(mxA);        
        mxArray *mxR = create_mat<T>(1, n);
        T *r = (T*)mxGetData(mxR);    
    
        for (int i = 0; i < n; ++i)
        {
            r[i] = det(A[i]);
        }
        
        return mxR;
    }
    else if (test_size(mxA, 4, 1, n))
    {        
        const Mat2x2<T>* A = (const Mat2x2<T>*)mxGetData(mxA);
        mxArray *mxR = create_mat<T>(1, n);
        T *r = (T*)mxGetData(mxR);   
        
        for (int i = 0; i < n; ++i)
        {
            r[i] = det(A[i]);
        }
        
        return mxR;
    }
    else if (test_size(mxA, 3, 1, n))
    {
        const SMat2x2<T>* A = (const SMat2x2<T>*)mxGetData(mxA);
        mxArray *mxR = create_mat<T>(1, n);
        T *r = (T*)mxGetData(mxR);
        
        for (int i = 0; i < n; ++i)
        {
            r[i] = det(A[i]);
        }
        
        return mxR;
    }
    else
    {
        mexErrMsgIdAndTxt("det2x2:invalidarg", "The size of input array is invalid.");
    }
    
    return NULL;
}



/**
 * Main entry
 *
 * Input:
 *   [0]: A  the input array containing the input matrices
 *           can be 2 x 2 x n, or 4 x n, or 3 x n.
 *
 * Output:
 *   [0]: R  the output vector (of size 1 x n) 
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgIdAndTxt("det2x2:invalidarg", 
                "The number of input arguments to det2x2 should be 1.");
    
    const mxArray *mxA = prhs[0];
    
    general_input_check(mxA, "det2x2:invalidarg");
        
    mxArray *mxR = 0;
    
    if (mxIsDouble(mxA))
    {
        mxR = do_det2x2<double>(mxA);
    }
    else if (mxIsSingle(mxA))
    {
        mxR = do_det2x2<float>(mxA);
    }
    else
    {
        mexErrMsgIdAndTxt("det2x2:invalidarg", 
                "The input array should be of class double or single.");
    }    
    
    plhs[0] = mxR;
}

