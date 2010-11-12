/********************************************************************
 *
 *  wmedian_cimp.cpp
 *
 *  The implementation of the core part of wmedian
 *
 *  Created by Dahua Lin, on Sep 28, 2010
 *
 ********************************************************************/

#include <mex.h>


template<typename T>
inline T find_median(int n, const T *x, const T *f)  // x: sorted values, // f: cumsum of sorted weights
{
    T t = f[n-1] / 2;
    
    int i = 0;
    while (i < n && f[i] < t) ++i;
    
    if (f[i] > t)
    {
        return x[i];
    }
    else
    {
        return i < n-1 ? (x[i] + x[i+1]) / 2 : x[i];
    }
}


template<typename T>
void find_medians(int m, int n, const T* X, const T *F, T *M)
{
    const T *x = X;
    const T *f = F;
    for (int i = 0; i < n; ++i, x += m, f += m)
    {
        M[i] = find_median(m, x, f);                
    }
}




// main entry
// Input
// [0]:     X   (sorted values) (along dim 1)
// [1]:     F   (cumsum of sorted weights)
// Output
// [0]:     M   (weighted median values)
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *mxX = prhs[0];
    const mxArray *mxF = prhs[1];
    
    int m = mxGetM(mxX);
    int n = mxGetN(mxX);
    
    mxArray *mxM = 0;
    
    if (mxIsDouble(mxX))
    {
        mxM = mxCreateNumericMatrix(1, n, mxDOUBLE_CLASS, mxREAL);
        
        find_medians(m, n, 
                (const double*)mxGetData(mxX), 
                (const double*)mxGetData(mxF), 
                (double*)mxGetData(mxM));
    }
    else // single
    {
        mxM = mxCreateNumericMatrix(1, n, mxSINGLE_CLASS, mxREAL);
        
        find_medians(m, n, 
                (const float*)mxGetData(mxX), 
                (const float*)mxGetData(mxF), 
                (float*)mxGetData(mxM));
    }
    
    plhs[0] = mxM;
}


