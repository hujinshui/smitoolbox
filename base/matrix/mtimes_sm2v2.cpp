/********************************************************************
 *
 *  mtimes_sm2v2.cpp
 *
 *  The C++ mex function for multiplication between a 2 x 2 symmetric
 *  matrix and a generic 2 x n matrix.
 *
 *  Created by Dahua Lin, on June 18, 2010
 *
 ********************************************************************/

#include "smallmat.h"


template<typename T>
mxArray* do_sm2v2(const mxArray *mxA, const mxArray *mxX, mxClassID cid)
{
    const SMat2x2<T> *A = (const SMat2x2<T>*)mxGetData(mxA);
    const Vec2<T>* X = (const Vec2<T>*)mxGetData(mxX);
    
    int m = (int)mxGetN(mxA);
    int n = (int)mxGetN(mxX);        
    
    mxArray *mxY = 0;
    
    if (m == 1)
    {
        mxY = mxCreateNumericMatrix(2, n, cid, mxREAL);
        const SMat2x2<T>& cA = A[0];
        Vec2<T>* Y = (Vec2<T>*)mxGetData(mxY);
        
        if (n == 1)
        {
            mtimes(cA, X[0], Y[0]);
        }
        else
        {
            for (int i = 0; i < n; ++i)
            {
                mtimes(cA, X[i], Y[i]);
            }
        }
    }
    else if (m == n)
    {
        mxY = mxCreateNumericMatrix(2, n, cid, mxREAL);
        Vec2<T>* Y = (Vec2<T>*)mxGetData(mxY);
        
        for (int i = 0; i < n; ++i)
        {
            mtimes(A[i], X[i], Y[i]);
        }
    }
    else
    {
        mexErrMsgIdAndTxt("mtimes_sm2v2:invalidarg", 
                "The sizes of A and X are inconsistent.");
    }
        
    return mxY;
}




/***
 * Main entry:
 *
 * Input
 *  [0]: A:     the 2x2 symmetric matrices [3x1 or 3xn]
 *  [1]: X:     the vectors [2xn]
 * Ouput
 *  [2]: Y:     the results
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxA = prhs[0];
    const mxArray *mxX = prhs[1];
    
    if (!(mxGetNumberOfDimensions(mxA) == 2 && !mxIsSparse(mxA) && mxGetM(mxA) == 3))
        mexErrMsgIdAndTxt("mtimes_sm2v2:invalidarg", 
                "A should be a non-sparse matrix with size(A,1) == 3.");
    
    if (!(mxGetNumberOfDimensions(mxX) == 2 && !mxIsSparse(mxX) && mxGetM(mxX) == 2))
        mexErrMsgIdAndTxt("mtimes_sm2v2:invalidarg", 
                "X should be a non-sparse matrix with size(X,1) == 2.");
        
    mxArray *mxY = 0;
    
    if (mxIsDouble(mxA) && mxIsDouble(mxX))
    {
        mxY = do_sm2v2<double>(mxA, mxX, mxDOUBLE_CLASS);
    }                
    else if (mxIsSingle(mxA) && mxIsSingle(mxX))
    {        
        mxY = do_sm2v2<float>(mxA, mxX, mxSINGLE_CLASS);
    }
    else
        mexErrMsgIdAndTxt("mtimes_sm2v2:invalidarg", 
                "A and X should be both of class double or both of class single.");
            
    plhs[0] = mxY;
}





