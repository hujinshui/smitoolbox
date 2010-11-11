/********************************************************************
 *
 *  intcount.cpp
 *
 *  The C++ mex implementation of integer counting
 *
 *  Created by Dahua Lin, on May 26, 2010
 *
 ********************************************************************/

#include <mex.h>

template<typename T>
void get_vs(const T *v, int ne, int& v0, int& v1)
{
    if (ne == 1)
    {
        v0 = 1;
        v1 = (int)(*v);
    }
    else
    {
        v0 = (int)(v[0]);
        v1 = (int)(v[1]);
    }
}



inline void get_range(const mxArray *mxRgn, int& v0, int& v1)
{
    int ne = mxGetNumberOfElements(mxRgn);
    
    if (!( (ne == 1 || ne == 2) && !mxIsSparse(mxRgn)))
    {
        mexErrMsgIdAndTxt("intcount:invalidarg", 
                "The range [v0 v1] should be a (non-sparse) pair.");
    }
    
    if (mxIsDouble(mxRgn))
    {
        const double *v = (const double*)mxGetData(mxRgn);
        get_vs(v, ne, v0, v1);
    }
    else if (mxIsSingle(mxRgn))
    {
        const float *v = (const float*)mxGetData(mxRgn);
        get_vs(v, ne, v0, v1);
    }
    else if (mxIsInt32(mxRgn))
    {
        const int *v = (const int*)mxGetData(mxRgn);
        get_vs(v, ne, v0, v1);
    }
    else
    {
        mexErrMsgIdAndTxt("intcount:invalidarg", 
                "The range [v0 v1] should be of class double, single, or int32");
    }
}


template<typename T>
void count(int v0, int v1, const T *v, int n, double *c)
{
    for (int i = 0; i < n; ++i)
    {
        int cv = (int)(v[i]);
        if (cv >= v0 && cv <= v1)
        {
            ++ c[cv - v0];
        }
    }
}


inline mxArray* do_count(int v0, int v1, const mxArray *mxVals)
{
    int m = v1 - v0 + 1;
    mxArray *mxCount = mxCreateDoubleMatrix(1, m, mxREAL);
    double *c = mxGetPr(mxCount);
    
    int n = mxGetNumberOfElements(mxVals);
    
    if (mxIsDouble(mxVals))
    {
        const double *v = (const double*)mxGetData(mxVals);
        count(v0, v1, v, n, c);
    }
    else if (mxIsSingle(mxVals))
    {
        const float *v = (const float*)mxGetData(mxVals);
        count(v0, v1, v, n, c);
    }
    else if (mxIsInt32(mxVals))
    {
        const int *v = (const int*)mxGetData(mxVals);
        count(v0, v1, v, n, c);
    }
    else
    {
        mexErrMsgIdAndTxt("intcount:invalidarg", 
                "The class of values should be either double, single, int32");
    }
    
    return mxCount;
}


/***********
 *
 *  Main entry
 *   
 *  Input:
 *    [0]:  the value range in form of [v0, v1] [pair (double|single|int32)]
 *    [1]:  the array of values to count [array (double|single|int32)]
 *
 *  Output:
 *    [0]:  the vector of counts (of size 1 x (v1 - v0 + 1)) [double array]
 *
 *  No 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    if (nrhs != 2)
        mexErrMsgIdAndTxt("intcount:invalidarg", 
                "The number of inputs to intcount should be 2.");
    
    const mxArray *mxRgn = prhs[0];
    const mxArray *mxVals = prhs[1];
    
    int v0, v1;
    get_range(mxRgn, v0, v1);    
    
    plhs[0] = do_count(v0, v1, mxVals);
}


