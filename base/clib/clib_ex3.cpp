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

#include "marray.h"
#include "sorting.h"

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
    Array<double> val(n);
    Array<int> idx(n);
    
    quicksort_asc(n, src, val.data(), idx.data());  
    
//     mexPrintf("Sorted result:\n");
//     mexPrintf("values: ");
//     printd(n, val.data(), true);
//     mexPrintf("indices: ");
//     printidx(n, idx.data(), true);
    
    if (is_sorted_asc(val.data(), n))
    {
        mexPrintf("The results is sorted.\n");        
    }
    else
    {
        mexPrintf("The results is NOT sorted !!!\n");
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



