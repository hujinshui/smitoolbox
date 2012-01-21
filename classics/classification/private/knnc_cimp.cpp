/********************************************************************
 *
 *  knnc_cimp.cpp
 *
 *  C++ mex implementation of KNNC
 *
 *  Created by Dahua Lin, on Jan 21, 2012
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

using namespace bcs;
using namespace bcs::matlab;

void perclass_accum_num(int M, int K, int n, const int32_t *L, double *A)
{
    for (int i = 0; i < n; ++i, L += K, A += M)
    {
        for (int k = 0; k < K; ++k)
        {
            int c = L[k];
            if (c >= 0 && c < K) A[c] += 1;
        }
    }
}


void perclass_accum_weight(int M, int K, int n, 
        const int32_t *L, const double *W, double *A)
{
    for (int i = 0; i < n; ++i, L += K, W += K, A += M)
    {
        for (int k = 0; k < K; ++k)
        {
            int c = L[k];
            if (c >= 0 && c < K) A[c] += W[k];
        }
    }
}



/**
 * main entry:
 *
 * Inputs:
 *  [0]:  M:    The number of classes [double scalar]
 *  [0]:  L:    selected neighbor labels [int32 K x n]
 *  [1]:  W:    selected neighbor weights [empty or double K x n]
 *
 * Outputs:
 *  [0]:  A:    the results [M x n]
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mM(prhs[0]);
    const_marray mL(prhs[1]);
    const_marray mW(prhs[2]);
    
    int M = (int)mM.get_scalar<double>();
    int K = (int)mL.nrows();
    int n = (int)mL.ncolumns();
    
    // prepare output
    
    marray mA = create_marray<double>(M, n);
        
    if (mW.isempty())
    {
        perclass_accum_num(M, K, n, mL.data<int32_t>(), 
                mA.data<double>());
    }
    else
    {
        perclass_accum_weight(M, K, n, 
                mL.data<int32_t>(), mW.data<int32_t>(), 
                mA.data<double>());
    }    
    
}

BCSMEX_MAINDEF
        

