/**********************************************************
 *
 *  chain_bp_cimp.cpp
 *
 *  The C++ mex implementation of chain_bp
 *
 *  Created by Dahua Lin, on Feb 1, 2012
 *
 **********************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <cmath>

using namespace bcs;
using namespace bcs::matlab;


struct FullSO
{
    FullSO(int K, const double *B)
    : _K(K), _B(B)
    {   
    }
    
    double operator() (int i, int j) const
    {
        return _B[i + j * _K];
    }                
    
    int _K;
    const double *_B;
};


struct SimpleSO
{
    SimpleSO(const double *B)
    : b0(B[0]), b1(B[1]) 
    {
    }
            
    double operator() (int i, int j) const
    {
        return i == j ? b0 : b1;
    }
        
    double b0;
    double b1;
};


inline void normalize_exp(int K, double *x)
{        
    double max_x = x[0];
    for (int i = 1; i < K; ++i) 
    {
        if (max_x < x[i]) max_x = x[i];
    }
    
    double s = 0;
    for (int i = 0; i < K; ++i)
    {
        s += (x[i] = std::exp(x[i] - max_x));
    }
    
    double c = 1.0 / s;
    for (int i = 0; i < K; ++i)
    {
        x[i] *= c;
    }
}

inline void normalize(int len, double *x)
{
    double s = 0;
    for (int i = 0; i < len; ++i) s += x[i];
    double c = 1.0 / s;
    for (int i = 0; i < len; ++i) x[i] *= c;
}


inline void update_u(int K, const double *M, const double *L, double *U)
{
    if (M)
    {
        for (int k = 0; k < K; ++k) U[k] = std::log(M[k]) + L[k];
    }
    else
    {
        for (int k = 0; k < K; ++k) U[k] = L[k];
    }
    
    double max_u = U[0];
    for (int k = 1; k < K; ++k)
    {
        if (max_u < U[k]) max_u = U[k];
    }
    
    for (int k = 0; k < K; ++k)
    {
        U[k] = std::exp(U[k] - max_u);
    }
}



template<class SecondOrdMat>
void chain_forward(int K, int n, 
        const double *LA, const SecondOrdMat& B, double *fmsg, double *U)
{
    // initialize
    
    update_u(K, 0, LA, U);
    
    // main loop
    for (int i = 1; i < n; ++i)
    {
        // compute fmsg
        
        for (int k = 0; k < K; ++k)
        {
            double v = 0;
            for (int l = 0; l < K; ++l)
            {
                v += U[l] * B(l, k);
            }
            fmsg[k] = v;
        }
        
        normalize(K, fmsg);
        
        // switch to next
        
        if (i < n-1)
        {
            LA += K;
            update_u(K, fmsg, LA, U);
            
            fmsg += K;
        }        
    }    
}


template<class SecondOrdMat>
void chain_backward(int K, int n, 
        const double *LA, const SecondOrdMat& B, double *bmsg, double *U)
{
    // initialize
    
    LA += (n - 1) * K;
    bmsg += (n - 2) * K;
    
    update_u(K, 0, LA, U);
    
    // main loop
    for (int i = n-1; i > 0; --i)
    {
        // compute bmsg 
        
        for (int l = 0; l < K; ++l)
        {
            for (int k = 0; k < K; ++k)
            {
                bmsg[k] += U[l] * B(k, l);
            }
        }
        
        normalize(K, bmsg);
        
        // switch to next
        
        if (i > 1)
        {
            LA -= K;
            update_u(K, bmsg, LA, U);
            
            bmsg -= K;
        }        
    }
}


void evaluate_belief(int K, int n, const double *LA, 
        const double *fmsg, const double *bmsg, double *mu)
{    
    // first one
    
    for (int k = 0; k < K; ++k)
    {
        mu[k] = LA[k] + std::log(bmsg[k]);
    }
    normalize_exp(K, mu);
    
    mu += K;
    LA += K;
    bmsg += K;
    
    // middle ones
    
    for (int i = 1; i < n-1; ++i)
    {
        for (int k = 0; k < K; ++k)
        {
            mu[k] = LA[k] + std::log(fmsg[k]) + std::log(bmsg[k]);
        }        
        normalize_exp(K, mu);
        
        mu += K;
        LA += K;       
        fmsg += K;
        bmsg += K;
    }
    
    // last one
    
    for (int k = 0; k < K; ++k)
    {
        mu[k] = LA[k] + std::log(fmsg[k]);
    }
    normalize_exp(K, mu);

}




template<class SecondOrdMat>
inline void do_chain_bp(int K, int n, 
        const double *LA, const SecondOrdMat& B, 
        double *mu, double *fmsg, double *bmsg)
{    
    chain_forward(K, n, LA, B, fmsg, mu);        
    chain_backward(K, n, LA, B, bmsg, mu);    
    evaluate_belief(K, n, LA, fmsg, bmsg, mu);
}



/**
 * main entry:
 *
 * Input
 *   [0]: LA:   the first-order log-potential [K x n]
 *   [1]: B:    the second-order potential [K x K or 1 x 2]
 *
 * Ouput
 *   [0]: mu:       the marginal distributions [K x n]
 *   [1]: fmsg:     the forward messages [K x (n-1)]
 *   [2]: bmsg: the backward messages [K x (n-1)]
 */
void bcsmex_main(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // take inputs
    
    const_marray mLA(prhs[0]);
    const_marray mB(prhs[1]);
    
    int K = (int)mLA.nrows();
    int n = (int)mLA.ncolumns();
    
    const double *LA = mLA.data<double>();
    
    // prepare output
    
    marray mMu = create_marray<double>(K, n);
    marray mFmsg = create_marray<double>(K, n-1);
    marray mBmsg = create_marray<double>(K, n-1);
    
    double *mu = mMu.data<double>();
    double *fmsg = mFmsg.data<double>();
    double *bmsg = mBmsg.data<double>();
    
    
    // main
    
    if (mB.nelems() == 2)
    {
        SimpleSO B(mB.data<double>());        
        do_chain_bp(K, n, LA, B, mu, fmsg, bmsg);
    }
    else
    {
        FullSO B(K, mB.data<double>());
        do_chain_bp(K, n, LA, B, mu, fmsg, bmsg);
    }
    
    // output
    
    plhs[0] = mMu.mx_ptr();
    plhs[1] = mFmsg.mx_ptr();
    plhs[2] = mBmsg.mx_ptr();    
    
}


BCSMEX_MAINDEF
        