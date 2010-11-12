/********************************************************************
 *
 *  ddsample_cimp.cpp
 *
 *  The C++ mex implementation for ddsample: sampling from
 *  discrete distributions
 *
 *  Created by Dahua Lin, on Nov 7, 2010
 *
 ********************************************************************/

#include <mex.h>
#include "ddsample.h"

using namespace msp;

/**
 *  K:  the number of classes
 *  m:  the number of different probability distributions
 *  n:  the number of samples to be drawn from each distribution
 *  F:  the cumulative distribution matrix [K x m column-major]
 *  V:  the uniformly distributed random values [n x m column-major]
 *  
 *  samples: the obtained samples [n x m column-major]
 */
void do_ddsample(int K, int m, int n, const double *F, const double *V, int *samples)
{
    if (m == 1)
    {
        if (n == 1)     // m == 1 & n == 1
        {            
            *samples = locate_interval(K, F, *V);
        }
        else            // m == 1 & n > 1
        {
            for (int i = 0; i < n; ++i)
            {
                samples[i] = locate_interval(K, F, V[i]);
            }
        }
    }
    else
    {
        if (n == 1)     // m > 1 & n == 1
        {
            for (int j = 0; j < m; ++j, F += K)
            {
                samples[j] = locate_interval(K, F, V[j]);
            }
        }
        else            // m > 1 & n > 1
        {
            for (int j = 0; j < m; ++j, F += K)
            {
                for (int i = 0; i < n; ++i)
                {
                    *(samples++) = locate_interval(K, F, *(V++));
                }
            }
        }
    }
}


mxArray *samples_to_matlab(int n, int m, int *samples)
{
    int N = n * m;
    mxArray *mx = mxCreateDoubleMatrix(n, m, mxREAL);
    double *r = mxGetPr(mx);
    
    for (int i = 0; i < N; ++i)
        r[i] = samples[i] + 1;
    
    return mx;    
}




/**
 * main entry:
 *
 * Input:
 *   [0]: F:    the cumulative distribution function [K x m double]
 *   [1]: V:    the uniformly distributed random variables [n x m double]
 * 
 * Output:
 *   [0]: s:    the samples [n x m double one-based]
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxF = prhs[0];
    const mxArray *mxV = prhs[1];
    
    double *F = mxGetPr(mxF);
    double *V = mxGetPr(mxV);
            
    int K = mxGetM(mxF);
    int m = mxGetN(mxF);
    int n = mxGetM(mxV);
    
    // do sampling
    
    int *s = new int[m * n];
    
    do_ddsample(K, m, n, F, V, s);
    
    // output
    
    plhs[0] = samples_to_matlab(n, m, s);
    
        
    // release memory
    
    delete[] s;
    
}

