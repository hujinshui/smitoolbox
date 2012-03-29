/********************************************************************
 *
 *  hmm_backward_cimp.cpp
 *
 *  The C++ mex implementation of hmm_backward
 *
 *  Created by Dahua Lin, on Feb 1, 2012
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <cmath>

using namespace bcs;
using namespace bcs::matlab;


inline void mult_mat_vec(int K, const double *T, const double *x, double *y)
{
    for (int i = 0; i < K; ++i) y[i] = 0;
    
    for (int j = 0; j < K; ++j)
    {
        double xj = x[j];
        
        for (int i = 0; i < K; ++i)
        {
            y[i] += *(T++) * xj;
        }
    }
}

inline double vmax(int len, const double *v)
{
    double m = v[0];    
    for (int i = 1; i < len; ++i) 
    {
        if (m < v[i]) m = v[i];
    }
    return m;
}



// core function

void do_backward(int K, int n, const double *T, const double *L, 
        const double *Lc, double *B)
{
    double *b = B + (n - 1) * K;
    
    for (int k = 0; k < K; ++k) b[k] = 1.0;
    
    const double *ll = L + (n - 1) * K;
    
    double *u = new double[K];
    
    for (int i = n-1; i > 0; --i)
    {
        double maxll = vmax(K, ll);
        
        for (int k = 0; k < K; ++k)
        {
            u[k] = b[k] * std::exp(ll[k] - maxll);
        }        
        b -= K;
        
        mult_mat_vec(K, T, u, b);
        
        double c = std::exp(maxll - Lc[i]);
        
        for (int k = 0; k < K; ++k)
        {
            b[k] *= c;
        }
        
        ll -= K;
    }
    
    delete[] u;    
}



/**
 * main entry:
 *
 * Input
 *   [0]: T:    the transition matrix 
 *   [1]: L:    the log-likelihood matrix
 *   [2]: Lc:   the log-conditional vector
 *
 * Ouput
 *   [0]: Beta:    The resultant matrix of beta values
 */
void bcsmex_main(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // take inputs
    
    const_marray mT(prhs[0]);
    const_marray mL(prhs[1]);
    const_marray mLc(prhs[2]);
    
    int K = (int)mT.nrows();
    int n = (int)mL.ncolumns();
    
    const double *T = mT.data<double>();
    const double *L = mL.data<double>();
    const double *Lc = mLc.data<double>();
    
    // prepare output
    
    marray mBeta = create_marray<double>(K, n);
    double *B = mBeta.data<double>();
    
    // main
    
    do_backward(K, n, T, L, Lc, B);
    
    // output
    
    plhs[0] = mBeta.mx_ptr();
    
}

BCSMEX_MAINDEF
        

