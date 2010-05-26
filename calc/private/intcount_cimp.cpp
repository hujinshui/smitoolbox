/********************************************************************
 *
 *  intcount_cimp.cpp
 *
 *  The C++ mex implementation of integer counting
 *
 *  Created by Dahua Lin, on May 26, 2010
 *
 ********************************************************************/

#include <mex.h>

void count(int v0, int v1, const int *v, int n, double *c)
{
    for (int i = 0; i < n; ++i)
    {
        if (v[i] >= v0 && v[i] <= v1)
        {
            ++ c[v[i] - v0];
        }
    }
}


/***********
 *
 *  Main entry
 *   
 *  Input:
 *    [0]:  the value range in form of [v0, v1] [double pair]
 *    [1]:  the array of values to count [int32 array]
 *
 *  Output:
 *    [0]:  the vector of counts (of size 1 x (v1 - v0 + 1)) [double array]
 *
 *  No 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxRgn = prhs[0];
    const mxArray *mxVals = prhs[1];
    
    const double *rgn = mxGetPr(mxRgn);
    int v0 = (int)(rgn[0]);
    int v1 = (int)(rgn[1]);
    
    const int *v = (const int*)mxGetData(mxVals);
    int n = mxGetNumberOfElements(mxVals);
    
    // main
    
    int m = v1 - v0 + 1;
    
    mxArray *mxCount = mxCreateDoubleMatrix(1, m, mxREAL);
    double *c = mxGetPr(mxCount);
    
    count(v0, v1, v, n, c);
    
    // output
    
    plhs[0] = mxCount;
}
