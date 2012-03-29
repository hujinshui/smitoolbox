/**********************************************************
 *
 *  chain_viterbi_cimp.cpp
 *
 *  The C++ mex implementation of chain_viterbi
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
    
    double operator() (int t, int i, int j) const
    {
        return _B[i + j * _K];
    }                
    
    int _K;
    const double *_B;
};


struct FullSO_W
{
    FullSO_W(int K, const double *B, const double *w)
    : _K(K), _B(B), _w(w)
    {   
    }
    
    double operator() (int t, int i, int j) const
    {
        return _B[i + j * _K] * _w[t];
    }                
    
    int _K;
    const double *_B;
    const double *_w;
};


struct SimpleSO
{
    SimpleSO(const double *B)
    : b0(B[0]), b1(B[1]) 
    {
    }
            
    double operator() (int t, int i, int j) const
    {
        return i == j ? b0 : b1;
    }
        
    double b0;
    double b1;
};


struct SimpleSO_W
{
    SimpleSO_W(const double *B, const double *w)
    : b0(B[0]), b1(B[1]), _w(w) 
    {
    }
            
    double operator() (int t, int i, int j) const
    {
        return (i == j ? b0 : b1) * _w[t];
    }
        
    double b0;
    double b1;
    const double *_w;
};



template<class SecondOrd>
void viterbi_forward(int K, int n, 
        const double *A0, const double *A, const SecondOrd& B, 
        double *V, int *R, double& final_v, int& final_s)
{
    // initialize
    
    if (A0)
    {
        for (int k = 0; k < K; ++k) V[k] = A0[k] + A[k];
    }
    else
    {
        for (int k = 0; k < K; ++k) V[k] = A[k];
    }
    
    // main loop
    
    double opt_v = 0;
    int opt_s = 0;
    
    for (int i = 1; i < n; ++i)
    {
        const double *U = V;
        V += K;
        A += K;
        
        for (int k = 0; k < K; ++k)
        {
            opt_v = U[0] + B(i-1, 0, k);
            opt_s = 0;
            
            for (int l = 1; l < K; ++l)
            {
                double v = U[l] + B(i-1, l, k);
                if (v > opt_v)
                {
                    opt_v = v;
                    opt_s = l;
                }
            }
            
            V[k] = opt_v + A[k];
            R[k] = opt_s;
        }
        
        R += K;
    }
    

    // finalize
    
    opt_v = V[0];
    opt_s = 0;
    
    for (int k = 1; k < K; ++k)
    {
        if (V[k] > opt_v)
        {
            opt_v = V[k];
            opt_s = k;
        }
    }
    
    final_v = opt_v;
    final_s = opt_s;
    
}


void viterbi_backtrace(int K, int n, int final_s, const int *R, 
        double *x)
{
    int s = final_s;
    x[n - 1] = s + 1;
    
    R += (n - 2) * K;
    
    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = ((s = R[s]) + 1);
        R -= K;
    }
}




template<class SecondOrd>
inline double do_viterbi(int K, int n, 
        const double *A0, const double *A, const SecondOrd& B, 
        double *x)
{
    double *V = new double[K * n];
    int *R = new int[K * (n-1)];
    
    double final_v;
    int final_s;
        
    viterbi_forward(K, n, A0, A, B, V, R, final_v, final_s);
    viterbi_backtrace(K, n, final_s, R, x);
    
    delete[] V;
    delete[] R;
    
    return final_v;
}



/**
 * main entry:
 *
 * Input
 *   [0]: A0:   additional first-order potential to the first node [empty or K x 1]
 *   [1]: A:    the first-order potential [K x n]
 *   [2]: B:    the second-order potential [K x K or 1 x 2]
 *   [3]: w:    the second-order link weights [empty or 1 x (n-1)]
 *
 * Ouput
 *   [0]: x:    the solved sequence
 *   [1]: v:    the optimal total potential
 */
void bcsmex_main(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // take inputs
    
    const_marray mA0(prhs[0]);
    const_marray mA(prhs[1]);
    const_marray mB(prhs[2]);
    const_marray mW(prhs[3]);
    
    int K = (int)mA.nrows();
    int n = (int)mA.ncolumns();
    
    const double *A0 = !mA0.is_empty() ? mA0.data<double>() : NULL;
    const double *A = mA.data<double>();
    const double *w = !mW.is_empty() ? mW.data<double>() : NULL;
    
    // prepare output
    
    marray mX = create_marray<double>(1, n);
    double *x = mX.data<double>();
    double v = 0;
    
    
    // main
    
    if (mB.nelems() == 2)
    {
        if (w)
        {
            SimpleSO_W B(mB.data<double>(), w);
            v = do_viterbi(K, n, A0, A, B, x);
        }
        else
        {
            SimpleSO B(mB.data<double>());
            v = do_viterbi(K, n, A0, A, B, x);
        }

    }
    else
    {
        if (w)
        {                        
            FullSO_W B(K, mB.data<double>(), w);
            v = do_viterbi(K, n, A0, A, B, x);
        }
        else
        {
            FullSO B(K, mB.data<double>());
            v = do_viterbi(K, n, A0, A, B, x);
        }
    }
    
    // output
    
    plhs[0] = mX.mx_ptr(); 
    plhs[1] = mxCreateDoubleScalar(v);
    
}


BCSMEX_MAINDEF






