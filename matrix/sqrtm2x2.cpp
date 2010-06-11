/********************************************************************
 *
 *  sqrtm2x2.cpp
 *
 *  The C++ mex implementation of matrix square root of 2x2 matrices
 *
 *  Created by Dahua Lin, on Jun 10, 2010
 *
 ********************************************************************/


#include "smallmat.h"


template<typename T>
inline mxArray* do_sqrtm2x2(const mxArray *mxA)
{
    int n = 0;
    
    if (test_size(mxA, 2, 2, n))
    {        
        const Mat2x2<T>* A = (const Mat2x2<T>*)mxGetData(mxA);
        mxArray *mxR = create_cube<T>(2, 2, n);
        Mat2x2<T>* R = (Mat2x2<T>*)mxGetData(mxR);
        
        for (int i = 0; i < n; ++i)
        {
            sqrtm(A[i], R[i]);
        }
        
        return mxR;
    }
    else if (test_size(mxA, 4, 1, n))
    {        
        const Mat2x2<T>* A = (const Mat2x2<T>*)mxGetData(mxA);
        mxArray *mxR = create_mat<T>(4, n);
        Mat2x2<T>* R = (Mat2x2<T>*)mxGetData(mxR);
        
        for (int i = 0; i < n; ++i)
        {
            sqrtm(A[i], R[i]);
        }
        
        return mxR;
    }
    else if (test_size(mxA, 3, 1, n))
    {
        const SMat2x2<T>* A = (const SMat2x2<T>*)mxGetData(mxA);
        mxArray *mxR = create_mat<T>(3, n);
        SMat2x2<T>* R = (SMat2x2<T>*)mxGetData(mxR);
        
        for (int i = 0; i < n; ++i)
        {
            sqrtm(A[i], R[i]);
        }
        
        return mxR;
    }
    else
    {
        mexErrMsgIdAndTxt("sqrtm2x2:invalidarg", "The size of input array is invalid.");
    }
}



/**
 * Main entry
 *
 * Input:
 *   [0]: A  the input array containing the input matrices
 *           can be 2 x 2 x n, or 4 x n, or 3 x n.
 *
 * Output:
 *   [0]: R  the output array (in the same size of A) 
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgIdAndTxt("sqrtm2x2:invalidarg", 
                "The number of input arguments to inv2x2 should be 1.");
    
    const mxArray *mxA = prhs[0];
    
    general_input_check(mxA, "sqrtm2x2:invalidarg");
        
    mxArray *mxR = 0;
    
    if (mxIsDouble(mxA))
    {
        mxR = do_sqrtm2x2<double>(mxA);
    }
    else if (mxIsSingle(mxA))
    {
        mxR = do_sqrtm2x2<float>(mxA);
    }
    else
    {
        mexErrMsgIdAndTxt("sqrtm2x2:invalidarg", 
                "The input array should be of class double or single.");
    }    
    
    plhs[0] = mxR;
}




