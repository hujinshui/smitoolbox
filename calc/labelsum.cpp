/********************************************************************
 *
 *  labelsum.cpp
 *
 *  The C++ mex implementation of labelsum.m
 *
 *  Created by Dahua Lin, on Jul 9, 2010
 *
 ********************************************************************/


#include <mex.h>


template<typename T>
void sum_columns(const T *X, T *Y, const int *L, int m, int n, int K)
{
    if (m == 1)
    {
        for (int j = 0; j < n; ++j)
        {
            Y[L[j]] += X[j];
        }
    }
    else
    {
        for (int j = 0; j < n; ++j)
        {            
            T *y = Y + L[j] * m;
            const T *x = X + j * m;
            
            for (int i = 0; i < m; ++i)
            {
                y[i] += x[i];
            }
        }
    }   
}


template<typename T>
void sum_rows(const T *X, T * Y, const int *L, int m, int n, int K)
{
    if (n == 1)
    {
        for (int i = 0; i < m; ++i)
        {
            Y[L[i]] += X[i];
        }
    }
    else
    {
        for (int i = 0; i < m; ++i)
        {
            T *y = Y + L[i];
            const T *x = X + i;
            
            for (int j = 0; j < n; ++j, x += m, y += K)
            {
                *y += *x;
            }
        }
    }
}


template<typename T>
inline void sum_vecs(const T *X, T *Y, const int *L, int m, int n, int K, bool is_sum_cols)
{
    if (is_sum_cols)
        sum_columns(X, Y, L, m, n, K);
    else
        sum_rows(X, Y, L, m, n, K);
}



template<typename T>
inline void copy_labels(const T *L0, int *L, int n, int K)
{
    for (int i = 0; i < n; ++i)
    {
        L[i] = (int)(L0[i]) - 1;  // from one-based to zero-based index
        
        if (L[i] < 0 || L[i] >= K)
            mexErrMsgIdAndTxt("labelsum:outofrange", 
                    "The given labels are out of range [1, K]");
    }
}




/*******
 * 
 * Input:
 *   [0]:  X:   the input matrix [m x n double|single|int32]
 *   [1]:  K:   the label range [double scalar]
 *   [2]:  L:   the label vector [m x 1 or 1 x n double|single|int32]
 *  
 * Output:
 *   [0]:  Y:   the output matrix [m x K or K x n double|single|int32]
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    if (nrhs != 3)
        mexErrMsgIdAndTxt("labelsum:invalidarg", 
                "The number of input arguments for labelsum should be 3.");
    
    const mxArray *mxX = prhs[0];
    const mxArray *mxK = prhs[1];
    const mxArray *mxL = prhs[2];
    
    // verify X
    
    if (!(mxGetNumberOfDimensions(mxX) == 2 && !mxIsSparse(mxX) && !mxIsComplex(mxX)))
        mexErrMsgIdAndTxt("labelsum:invalidarg", 
                "X should be a non-sparse real matrix.");
    
    int m = mxGetM(mxX);
    int n = mxGetN(mxX);
    
    if (!(mxGetNumberOfDimensions(mxL) == 2 && !mxIsSparse(mxL)))
        mexErrMsgIdAndTxt("labelsum:invalidarg", 
                "L should be a non-sparse matrix.");
    
    
    // verify K
    
    if (mxGetNumberOfElements(mxK) != 1 && !mxIsSparse(mxK))
        mexErrMsgIdAndTxt("labelsum:invalidarg", 
                "K should be a non-sparse scalar.");
    
    mxClassID K_cid = mxGetClassID(mxK);
    
    int K = 0;
    
    switch (K_cid)
    {
        case mxDOUBLE_CLASS:
            K = (int)(*(mxGetPr(mxK)));
            break;
        case mxSINGLE_CLASS:
            K = (int)(*((float*)mxGetData(mxK)));
            break;
        case mxINT32_CLASS:
            K = *((int*)mxGetData(mxK));
            break;
        default:
            mexErrMsgIdAndTxt("labelsum:invalidarg", 
                    "K should be of class double, single, or int32");            
    }
    
    if (K < 1)
        mexErrMsgIdAndTxt("labelsum:invalidarg", 
                "K should be a positive value not less than 1.");    
        
    // verify L
    
    int mL = mxGetM(mxL);
    int nL = mxGetN(mxL);
    
    bool is_sum_cols;
    int numL;
    
    if (mL == 1 && nL == n)
    {
        is_sum_cols = true;
        numL = n;
    }
    else if (mL == m && nL == 1)
    {
        is_sum_cols = false;
        numL = m;
    }
    else
    {
        mexErrMsgIdAndTxt("labelsum:invalidarg", 
                "The size of L should be either 1 x n or m x 1");
    }
        
    // copy labels
    
    int *L = new int[numL];        
    
    switch (mxGetClassID(mxL))
    {
        case mxDOUBLE_CLASS:
            copy_labels((double*)mxGetData(mxL), L, numL, K);
            break;
        case mxSINGLE_CLASS:
            copy_labels((float*)mxGetData(mxL), L, numL, K);
            break;
        case mxINT32_CLASS:
            copy_labels((int*)mxGetData(mxL), L, numL, K);
            break;
        default:
            mexErrMsgIdAndTxt("labelsum:invalidarg", 
                    "L should be of class double, single, or int32");
    }
    
    
    // do the computation
    
    mxArray *mxY = 0;
    
    int mY, nY;
    if (is_sum_cols)
    {
        mY = m;
        nY = K;
    }
    else
    {
        mY = K;
        nY = n;
    }
    
    switch (mxGetClassID(mxX))
    {
        case mxDOUBLE_CLASS:
            mxY = mxCreateDoubleMatrix(mY, nY, mxREAL);
            sum_vecs((const double*)mxGetData(mxX), (double*)mxGetData(mxY), 
                    L, m, n, K, is_sum_cols);
            break;
        case mxSINGLE_CLASS:
            mxY = mxCreateNumericMatrix(mY, nY, mxSINGLE_CLASS, mxREAL);
            sum_vecs((const float*)mxGetData(mxX), (float*)mxGetData(mxY), 
                    L, m, n, K, is_sum_cols);
            break;
        case mxINT32_CLASS:
            mxY = mxCreateNumericMatrix(mY, nY, mxINT32_CLASS, mxREAL);
            sum_vecs((const int*)mxGetData(mxX), (int*)mxGetData(mxY), 
                    L, m, n, K, is_sum_cols);
            break;
        default:
            mexErrMsgIdAndTxt("labelsum:invalidarg", 
                    "X should be of class double, single, or int32");  
    }   
    
    
    // release memory
    
    delete[] L;
    
    
    // output
    plhs[0] = mxY;
}





