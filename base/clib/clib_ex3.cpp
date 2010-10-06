/********************************************************************
 *
 *  clib_ex3.cpp
 *
 *  The example shows how to use sort functions
 *
 *  Created by Dahua Lin, on Oct 5, 2010
 *
 ********************************************************************/

// #define MATLAB_DUMP_QSORT_ACTION

#define _SECURE_SCL 0

#include "marray.h"
#include <vector>
#include <algorithm>

using namespace smi;

void printidx(int n, const int *v, bool ret = false)
{
    for (int i = 0; i < n; ++i) mexPrintf("%d ", v[i]+1);
    if (ret) mexPrintf("\n");
}

void printd(int n, const double *v, bool ret = false)
{
    for (int i = 0; i < n; ++i) mexPrintf("%.4g ", v[i]);
    if (ret) mexPrintf("\n");
}


void run_on(int n, const double *src)
{
    std::vector<double> v(n);
    
    std::copy(src, src+n, v.begin());
    
    std::sort(v.begin(), v.end());
    
    for (int i = 0; i < n; ++i)
    {
        mexPrintf("%.4g ", v[i]);
    }    
}


/**
 * Main entry:
 *
 * Input:
 *   [0]: src  the source array (should be a double vector)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgTxt("The number of inputs should be exactly 1.");
    
    MArray mSrc(prhs[0]);
    
    if (!(mSrc.is_double() && !mSrc.is_sparse()))
    {
        mexErrMsgTxt("The src should be a non-sparse double vector.");
    }
        
    run_on(mSrc.nelems(), mSrc.get_data<double>());
    mexPrintf("\n");
}



