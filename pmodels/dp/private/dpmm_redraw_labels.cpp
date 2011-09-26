/**********************************************************
 *
 *  dpmm_redraw_labels.cpp
 *
 *  C++ mex implementation of label re-drawing
 *
 *  Created by Dahua Lin, on Sep 24, 2011
 *
 **********************************************************/

#include <bcslib/matlab/bcs_mex.h>

#include "dp_base.h"
#include <cmath>

using namespace bcs;
using namespace bcs::matlab;

inline double draw_label(int K, const double *e, double e0, double u, double *w)
{
    double mv = e0;
    for (int k = 0; k < K; ++k)
    {
        if (e[k] > mv) mv = e[k];
    }
    
    double w0 = std::exp(e0 - mv);
    double sw = w0;
    
    for (int k = 0; k < K; ++k)
    {
        sw += (w[k] = std::exp(e[k] - mv));
    }
    
    int z = dd_draw(K, w, sw, u);
    return z < K ? double(z + 1) : 0.0;
}


/**
 * Input
 *   [0] E:     the main log-evidence matrix [K x n]
 *   [1] ev0:   the base log-evidence vector [1 x n]
 *   [2] rnums: the random numbers [1 x n]
 * 
 * Output
 *   [0] z:     the obtained labels [1 x n]
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_marray mE(prhs[0]);
    const_marray mEv0(prhs[1]);
    const_marray mRnums(prhs[2]);
    
    int K = (int)mE.nrows();
    int n = (int)mE.ncolumns();
    
    const double *E = mE.data<double>();
    const double *ev0 = mEv0.data<double>();
    const double *rnums = mRnums.data<double>();
    
    marray mZ = create_marray<double>(1, n);
    double *z = mZ.data<double>();
    
    // main
    
    double *w = new double[K];
    
    for (int i = 0; i < n; ++i)
    {
        z[i] = draw_label(K, E + K * i, ev0[i], rnums[i], w);
    }
    
    delete[] w;
        
    // output
    
    plhs[0] = mZ.mx_ptr();
}

BCSMEX_MAINDEF

