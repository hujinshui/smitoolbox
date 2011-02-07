/********************************************************************
 *
 *  gatrwa_cimp.cpp
 *
 *  The C++ mex implementation of gatrwa_cimp.m
 *
 *  Created by Dahua Lin, on Feb 6, 2011
 *
 ********************************************************************/


#include "../../../graph/clib/graph_mex.h"

#include <cmath>

using namespace smi;

typedef CRefAdjList<double, boost::undirected_tag> graph_t;


inline double calc_b(vertex_t v, const graph_t& eg, int m, 
        const double *sigma, const double *rho)
{
    double b = 0;
    graph_t::out_edge_iterator ep, eend;
    for (boost::tie(ep, eend) = out_edges(v, eg); ep != eend; ++ep)
    {
        edge_t e = *ep;
        double w = eg.get_weight(e);
        
        int ei = e.i < m ? e.i : e.i - m;
        double r = rho[ei];
        
        vertex_t vt = target(e, eg);
        double s = sigma[vt.i];
        
        b += w * r * s;
    }
    
    return b;
}


inline void update_sigma(double *sigma, vertex_t v, const graph_t& eg, int m, 
        const double *Jdv, const double *rho)
{
    double Jv = Jdv[v.i];    
    double b = calc_b(v, eg, m, sigma, rho);        
    double sig = (std::sqrt(b*b + 4 * Jv) - b) / (2 * Jv);
    
    sigma[v.i] = sig;
}



void do_compute_bvec(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mEg(prhs[1]);
    MArray mSigma(prhs[2]);
    MArray mRho(prhs[3]);
    
    matlab_graph_repr egr(mEg);    
    if (egr.weight_class() != mxDOUBLE_CLASS)
    {
        mexErrMsgIdAndTxt("gatrwa:invalidarg", 
                "The edge weights should be of double class.");
    }  
    graph_t eg = egr.to_cref_wadjlist_ud<double>();
    int n = (int)num_vertices(eg);
    int m = (int)num_edges(eg);
    
    const double *sigma = mSigma.get_data<double>();
    const double *rho = mRho.get_data<double>();
    
    // compute
    
    mxArray *mxB = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *b = mxGetPr(mxB);
    
    for (int i = 0; i < n; ++i)
    {
        b[i] = calc_b(i, eg, m, sigma, rho);
    }    
    
    plhs[0] = mxB;
    
}

void do_update_sigma(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mEg(prhs[1]);
    MArray mJdv(prhs[2]);
    mxArray *mxSigma = mxDuplicateArray(prhs[3]);
    MArray mRho(prhs[4]);
    MArray mOrd(prhs[5]);
    
    matlab_graph_repr egr(mEg);    
    if (egr.weight_class() != mxDOUBLE_CLASS)
    {
        mexErrMsgIdAndTxt("gatrwa:invalidarg", 
                "The edge weights should be of double class.");
    }  
    
    
    graph_t eg = egr.to_cref_wadjlist_ud<double>();
    
    const double *Jdv = mJdv.get_data<double>();
    double *sigma = mxGetPr(mxSigma);
    const double *rho = mRho.get_data<double>();
    const int *ord = mOrd.get_data<int>();
    
    int N = mOrd.nelems();
    int m = (int)num_edges(eg);
    
    // compute
    
    for (int i = 0; i < N; ++i)
    {
        int v = ord[i];
        
        update_sigma(sigma, v, eg, m, Jdv, rho);
    }
    
    plhs[0] = mxSigma;
}


/**
 * main entry:
 *
 * Input
 *   [0]: code:  the operation code [double scalar]
 *
 * If code == 0:  compute the b-vector
 *
 * Input
 *      [1]: eg:        obj.egraph
 *      [2]: sigma:     obj.sigma;
 *      [3]: rho:       obj.rho
 * Output
 *      [0]: b:         the b vector
 *
 *
 * If code == 1:  update sigma(s)
 *
 *  Input
 *      [1]: eg:        obj.egraph
 *      [2]: Jdv:       obj.Jdv
 *      [3]: sigma:     obj.sigma;
 *      [4]: rho:       obj.rho;
 *      [5]: ord:       order of updating
 *  Output:
 *      [0]: sigma:     updated sigma vector
 * 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    MArray mCode(prhs[0]);
    int code = (int)mCode.get_double_scalar();
    
    if (code == 0)
    {
        do_compute_bvec(nlhs, plhs, nrhs, prhs);
    }
    else if (code == 1)
    {
        do_update_sigma(nlhs, plhs, nrhs, prhs);
    }                    
}


