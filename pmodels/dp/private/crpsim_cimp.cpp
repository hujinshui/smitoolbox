/********************************************************************
 *
 *  crpsim_cimp.cpp
 *
 *  The C++ mex implementation of crpsim (CRP simulation)
 *
 *  Created by Dahua Lin, on Sep 17, 2011
 *
 ********************************************************************/


#include <bcslib/matlab/bcs_mex.h>
#include <vector>

#include "dp_base.h"

using namespace bcs;
using namespace bcs::matlab;

inline double calc_sum(int n, const double *a)
{
    double s(0);
    for (int i = 0; i < n; ++i) s += a[i];
    return s;
}


void do_sim(double alpha, int n, double *x, std::vector<double>& acc, 
        const double *randnums)
{
    int K = (int)acc.size();
    double tw = calc_sum(K, &(acc[0])) + alpha;
    
    for (int i = 0; i < n; ++i)
    {
        int k = dd_draw(K, &(acc[0]), tw, randnums[i]);
        
        if (k < K)
        {
            acc[k] += 1.0;
        }
        else
        {
            acc.push_back(1.0);
            ++ K;
        }       
        
        tw += 1.0;
        x[i] = double(k + 1);
    }
}



/***********
 *
 *  Main entry
 *   
 *  Input:
 *    [0] alpha:        the concentration parameter [double scalar]
 *    [1] randnums:     the sequence of uniform random numbers [1 x n]
 *    [2] pricount:     the prior count [double array]
 *
 *  Output:
 *    [0] x:            the generated sequence
 *    [1] acc:          the accumulated count of atoms
 *
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_marray mAlpha(prhs[0]);
    const_marray mRandNums(prhs[1]);
    const_marray mPriCount(prhs[2]);

    double alpha = mAlpha.get_scalar<double>();
    
    int n = (int)mRandNums.nelems();
    const double *randnums = mRandNums.data<double>();    
    
    int K = (int)mPriCount.nelems();  
    const double *pric = 0;
    if (K > 0)
    {
        pric = mPriCount.data<double>();
    }

    // prepare output
    
    marray mX = create_marray<double>(1, (size_t)n);
    double *x = mX.data<double>();
    
    // do simulation
    
    std::vector<double> acc;
    acc.reserve((size_t)(K + 16));
    
    if (K > 0)
    {
        for (int k = 0; k < K; ++k)
        {
            acc.push_back(pric[k]);
        }
    }
    
    do_sim(alpha, n, x, acc, randnums);
    
    // make output
    
    plhs[0] = mX.mx_ptr();   
    
    if (nlhs >= 2)
    {
        plhs[1] = to_matlab_row(acc).mx_ptr();
    }
}


BCSMEX_MAINDEF
        
        