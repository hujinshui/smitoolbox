/********************************************************************
 *
 *  dpmm_update_labels.h
 *
 *  Created by Dahua Lin, on Sep 20, 2011
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <cmath>

#include "dp_base.h"

using namespace bcs;
using namespace bcs::matlab;

inline int draw_label(
        int K, 
        double alpha, 
        double loglik0, 
        const double *logliks,      /* length = K */
        const double *pricounts,    /* length = K */
        const double *counts,       /* length = K */
        double rnum,
        double *temp                /* length = K + 1 */
        )
{
    // collect evidences in log-scale
    
    double t0 = loglik0 + std::log(alpha);
    double mv = t0;
    
    for (int k = 0; k < K; ++k)
    {
        double tk = logliks[k] + std::log(pricounts[k] + counts[k]);
        temp[k] = tk;
        if (mv > tk) mv = tk;        
    }
    temp[K] = t0;
    
    // shift and calculate sum-weight
    
    double sw = 0;
    for (int k = 0; k <= K; ++k)
    {
        sw += (temp[k] = std::exp(temp[k] - mv));
    }
    
    // draw sample
    
    return dd_draw(K, temp, sw, rnum);
    
}



/**
 * @param K         the number of active atoms
 * @param n         the number of observations
 * @param alpha     the concentration coefficient
 * @param logliks0  the 1 x n log marginal likelihood with respect to base
 * @param logliks   the K x n log-likelihood value matrix
 * @param pricounts the prior count vector of length K
 * @param counts    the current count vector of length K
 * @param labels    the label vector of length n (one-based)
 * @param len       the total length of sequence of indices to be updated
 * @param inds      the sequence of indices to be updated
 * @param cp        the current position (as input: starting, as output: next starting)
 */
bool do_main(int K, int n, double alpha, 
        const double *logliks0, const double *logliks, 
        const double *pricounts, double *counts, double *labels, 
        int len, const double *inds, const double *rnums, int& cp)
{
    double *temp = new double[K];
    
    bool ch = false;
    
    while (cp < len)
    {
        double rn = rnums[cp];
        int i = (int)inds[cp] - 1;
       
        // withdraw current label
        
        int z0 = int(labels[i] - 1);
        
        labels[i] = 0;
        if (z0 >= 0 && z0 < K)
        {
            -- counts[z0];
        }                
        
        // draw new labels
        
        int z = draw_label(K, alpha, logliks0[i], logliks + K * i,
                pricounts, counts, rn, temp);                
        
        if (z != z0) ch = true;
        
        if (z < K)
        {
            labels[i] = double(z+1);
            ++ counts[z];
            ++ cp;
        }
    }
    
    delete[] temp;
    
    return ch;
}



/**
 * main entry:
 *
 * Inputs
 *    [0] alpha:        the concentration parameter (double scalar)
 *    [1] logliks0:     the vector of log-margin-likelihood w.r.t base
 *    [2] logliks:      the matrix of log-likelihood (double K x n)
 *    [3] pricounts:    the prior counts of atoms (vector of len K)
 *
 * Inputs (can be updated)
 *    [4] counts:       the counts at current phases [vector of len K]
 *    [5] labels:       the labels [vector of len n]
 *
 * Inputs
 *    [6] inds:         the sequence of indices of the samples to be updated
 *    [7] rnums:        the sequence of random numbers 
 *    [8] sp:           the starting position (inclusive)
 *
 * Outputs
 *    [0] ep:           the ending position (inclusive)
 *    [1] ch:           whether changes to labels were made
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mAlpha(prhs[0]);
    const_marray mLogliks0(prhs[1]);
    const_marray mLogliks(prhs[2]);
    const_marray mPriCounts(prhs[3]);
    
    marray mCounts(const_cast<mxArray*>(prhs[4]));
    marray mLabels(const_cast<mxArray*>(prhs[5]));
    
    const_marray mInds(prhs[6]);
    const_marray mRNums(prhs[7]);
    const_marray mSp(prhs[8]);
    
    double alpha = mAlpha.get_scalar<double>();
    int K = (int)mLogliks.nrows();
    int n = (int)mLogliks.ncolumns();
    
    const double *logliks0 = mLogliks0.data<double>();
    const double *logliks = mLogliks.data<double>();
    const double *pricounts = mPriCounts.data<double>();
    
    double *counts = mCounts.data<double>();
    double *labels = mLabels.data<double>();
    
    int len = (int)mInds.nelems();
    const double *inds = mInds.data<double>();
    const double *rnums = mRNums.data<double>();
    int cp = (int)mSp.get_scalar<double>() - 1;
    
    // main
    
    int ep = 0;    
    bool ch = do_main(K, n, alpha, logliks0, logliks, pricounts, 
            counts, labels, len, inds, rnums, cp);
    
    // make outputs
    
    marray mEP = create_mscalar(double(cp + 1));
    marray mCh = create_mscalar(ch);
    
    plhs[0] = mEP.mx_ptr();
    plhs[1] = mCh.mx_ptr();
    
}


BCSMEX_MAINDEF
        
