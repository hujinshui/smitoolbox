/********************************************************************
 *
 *  repnum_cimp.cpp
 *
 *  The C++ mex implementation for repnum
 *
 *  Created by Dahua Lin, on June 6, 2010
 *
 ********************************************************************/

#include <mex.h>

template<typename T>
int get_nums(const mxArray *mxNs, int n, int *iv)
{
    const T *v = (const T*)mxGetData(mxNs);
    
    int s = 0;
    for (int i = 0; i < n; ++i)
    {
        s += (iv[i] = (int)(v[i]));
    }
    return s;
}


void do_repnum(int n, const int *iv, double *r)
{
    for (int k = 1; k <= n; ++k)
    {
        int c = iv[k-1];
        for (int i = 0; i < c; ++i) *r++ = k;
    }
}


/**
 *
 * main entry:
 *
 *  Input:
 *    [0]: ns:  the numbers of repeating times
 *
 *  Output:
 *    [1]: r:   the generated vector 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgIdAndTxt("repnum:invalidarg", 
                "The number of inputs to repnum should be 1.");
    
    const mxArray *mxNs = prhs[0];
        
    
    int ndims_ns = mxGetNumberOfDimensions(mxNs);
    int m_ns = mxGetM(mxNs);
    int n_ns = mxGetN(mxNs);
    
    if (!(!mxIsSparse(mxNs) && !mxIsComplex(mxNs) && !mxIsEmpty(mxNs) && 
            ndims_ns == 2 && (m_ns == 1 || n_ns == 1)))
        mexErrMsgIdAndTxt("repnum:invalidarg", 
                "ns should be a non-empty real full vector.");
    
    int n = m_ns * n_ns;
    int *iv = new int[n];
    int N = 0;
    
    switch (mxGetClassID(mxNs))
    {
        case mxDOUBLE_CLASS:
            N = get_nums<double>(mxNs, n, iv);
            break;
        case mxSINGLE_CLASS:
            N = get_nums<float>(mxNs, n, iv);
            break;
        case mxINT32_CLASS:
            N = get_nums<int>(mxNs, n, iv);
            break;
        default:
            mexErrMsgIdAndTxt("repnum:invalidarg", 
                    "ns should be of either of the class double, single, or int32.");
                
    }
    
    mxArray *mxR = m_ns == 1 ? 
        mxCreateDoubleMatrix(1, N, mxREAL) :    
        mxCreateDoubleMatrix(N, 1, mxREAL);
        
    double *r = mxGetPr(mxR);    
    do_repnum(n, iv, r);
                    
    delete[] iv;
    
    // output
    
    plhs[0] = mxR;
        
}

