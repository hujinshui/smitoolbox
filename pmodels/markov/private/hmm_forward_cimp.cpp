/********************************************************************
 *
 *  hmm_forward_cimp.cpp
 *
 *  The C++ mex implementation of hmm_forward
 *
 *  Created by Dahua Lin, on Feb 1, 2012
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <cmath>

using namespace bcs;
using namespace bcs::matlab;


inline void nrm_process(int K, const double *P0, const double *L, 
        double *A, double& Lc)
{ 
    double maxL = L[0];
    for (int k = 1; k < K; ++k)
    {
        if (maxL < L[k]) maxL = L[k];
    }
    
    double s = 0;
    for (int k = 0; k < K; ++k)
    {
        s += (A[k] = P0[k] * std::exp(L[k] - maxL));
    }
    
    double c = 1.0 / s;    
    for (int k = 0; k < K; ++k)
    {
        A[k] *= c;
    }
    
    Lc = maxL + std::log(s);    
}


// core function

void do_forward(int K, int n, 
        const double *Pi, const double *T, const double *L, 
        double *A, double *Lc)
{
    // initialize
        
    nrm_process(K, Pi, L, A, Lc[0]);
        
    const double *A0 = A;
    A += K;
    L += K;
        
    // recursively proceed
    
    for (int i = 1; i < n; ++i)
    {        
        for (int k = 0; k < K; ++k)
        {
            double v = 0;            
            for (int u = 0; u < K; ++u)
            {
                v += A0[u] * T[u + k * K];
            }
            A[k] = v;
        }
        
        nrm_process(K, A, L, A, Lc[i]);
                
        A0 = A;
        A += K;
        L += K;
    }
     
}



/**
 * main entry:
 *
 * Input
 *   [0]: Pi:   the initial distribution
 *   [1]: T:    the transition matrix 
 *   [2]: L:    the log-likelihood matrix
 *
 * Ouput
 *   [0]: Alpha:    The resultant matrix of alpha values
 */
void bcsmex_main(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // take inputs
    
    const_marray mPi(prhs[0]);
    const_marray mT(prhs[1]);
    const_marray mL(prhs[2]);
    
    int K = (int)mT.nrows();
    int n = (int)mL.ncolumns();
    
    const double *Pi = mPi.data<double>();
    const double *T = mT.data<double>();
    const double *L = mL.data<double>();
    
    // prepare output
    
    marray mAlpha = create_marray<double>(K, n);
    double *A = mAlpha.data<double>();
    
    marray mLc = create_marray<double>(1, n);
    double *Lc = mLc.data<double>();
    
    // main
    
    do_forward(K, n, Pi, T, L, A, Lc);
    
    // output
    
    plhs[0] = mAlpha.mx_ptr();
    plhs[1] = mLc.mx_ptr();
    
}

BCSMEX_MAINDEF
        

