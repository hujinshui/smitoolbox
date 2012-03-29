/**********************************************************
 *
 *  topic_loglik_cimp.cpp
 *
 *  C++ mex implementation of topic_loglik
 *
 *  Created by Dahua Lin, on Feb 18, 2012
 *
 **********************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include "specfuncs.h"
#include "corpus.h"

using namespace bcs;
using namespace bcs::matlab;
using namespace smi;


double calc_loglik(int V, int K, const double *Beta, 
        const Doc& doc, const double *q, double *temp)
{
    int nw = doc.nwords;
    
    double ll = 0;
    
    for (int i = 0; i < nw; ++i)  // for each word
    {
        // calculate per-word posterior of topics
        
        int32_t wd = doc.words[i];
        double c = doc.counts[i];
        
        double s = 0;
        for (int k = 0; k < K; ++k) 
        {
            s += (temp[k] = q[k] * Beta[wd + V * k]);
        }
        
        double inv_s = 1.0 / s;
        for (int k = 0; k < K; ++k) temp[k] *= inv_s;
        
        // calculate generative probability of this word
        // with z marginalized out
        
        double pv = 0;
        for (int k = 0; k < K; ++k)
        {
            pv += temp[k] * Beta[wd + V * k];
        }
        
        // add to total loglik
        
        ll += c * std::log(pv);        
    }
    
    return ll;
}



/**
 * Inputs
 *   [0] Beta:          The word distributions of topics [V x K]
 *   [1] I:             The whole list of words [int32 zero-based]
 *   [2] J:             The whole list of documents [int32 zero-based]
 *   [3] C:             The whole list of word counts [double]
 *   [4] Q:             The per-doc posterior distr. of topics [K x n]
 *
 * Outputs
 *   [0] v:             The per-doc log-likelihood values [1 x n]
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mBeta(prhs[0]);
    const_marray mI(prhs[1]);
    const_marray mJ(prhs[2]);
    const_marray mC(prhs[3]);
    const_marray mQ(prhs[4]);
    
    int V = mBeta.nrows();
    int K = mBeta.ncolumns();
    int n = mQ.ncolumns();
    int len = mI.nelems();
    
    const double *Beta = mBeta.data<double>();
    const int32_t *I = mI.data<int32_t>();
    const int32_t *J = mJ.data<int32_t>();
    const double *C = mC.data<double>();
    const double *Q = mQ.data<double>();
    
    // prepare output
    
    marray mV = create_marray<double>(1, n);
    double *v = mV.data<double>();
    
    // main
    
    Corpus corpus(n, len, I, J, C);
    const double *q = Q;
    double *temp = new double[K];
    
    for (int i = 0; i < n; ++i, q += K)
    {
        const Doc& doc = corpus.doc(i);
        v[i] = calc_loglik(V, K, Beta, doc, q, temp);
    }
    
    delete[] temp;
    
    // output
    
    plhs[0] = mV.mx_ptr();
        
}


BCSMEX_MAINDEF
        
        