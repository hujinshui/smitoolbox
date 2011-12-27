/********************************************************************
 *
 *  wishartd_sample_cimp.cpp
 *
 *  The C++ mex implementation for pieces
 *
 *  Created by Dahua Lin, on Mar 22, 2011
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <cmath>

using namespace bcs;
using namespace bcs::matlab;

// core function

void wishartd_gen(int d, const double *css, const double *nrms, double *R)
{
    double *dv = R;
    for (int i = 0; i < d; ++i, dv += (d+1))
    {
        *dv = std::sqrt(*(css++));
    }
    
    for (int i = 1; i < d; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            R[i + j * d] = *(nrms++);
        }
    }
}


/**
 * The main entry
 *
 *  Input:
 *      [1] css:        the chi-squared random numbers [d x n double]
 *      [2] nrms:       the standard normal random numbers [d(d-1)/2 x n double] 
 * 
 *  Output:
 *      [0] R:          the output [d x d x n]
 *
 *  Note:
 *      when d == 2 or d == 3, R store the sample matrix, 
 *
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_marray mCss(prhs[0]);
    const_marray mNrms(prhs[1]);           
    
    int d = (int)mCss.nrows();
    int n = (int)mCss.ncolumns();
    
    const double *css = mCss.data<double>();
    const double *nrms = mNrms.data<double>();
    
    // prepare output
    
    mxArray *mxR = 0;
    if (n == 1)
    {
        mxR = mxCreateDoubleMatrix((mwSize)d, (mwSize)d, mxREAL);
    }
    else
    {
        mwSize dims[3];
        dims[0] = (mwSize)d;
        dims[1] = (mwSize)d;
        dims[2] = (mwSize)n;
        
        mxR = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    }
    
    double *R = mxGetPr(mxR);
                    
    // main
        
    if (n == 1)
    {
        wishartd_gen(d, css, nrms, R);
    }
    else
    {
        int dn = d * (d - 1) / 2;
        int dr = d * d;
        
        for (int i = 0; i < n; ++i)
        {                        
            wishartd_gen(d, css + d * i, nrms + dn * i, R + dr * i);
        }
    }
       
    // output
    
    plhs[0] = mxR;
    
}


BCSMEX_MAINDEF
        
        