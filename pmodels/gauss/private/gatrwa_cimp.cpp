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


inline double compute_sigma(const graph_t& eg, int m, vertex_t v, 
        const double *Jdv, const double *sigma, const double *rho)
{
    double Jv = Jdv[v.i];    
    double b = calc_b(v, eg, m, sigma, rho);        
    double sig = (std::sqrt(b*b + 4 * Jv) - b) / (2 * Jv);
    
    return sig;
}

inline double compute_rho(const graph_t& eg, edge_t e, 
        const double *sigma, const double *ep)
{
    vertex_t s = source(e, eg);
    vertex_t t = target(e, eg);
    
    double w = eg.get_weight(e);
    double a = w * sigma[s.i] * sigma[t.i];
    
    double beta = ep[e.i];
    
    double r = (beta - std::sqrt(beta*beta + 4 * a*a)) / (2 * a);
    return r;
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
        sigma[i] = compute_sigma(eg, m, v, Jdv, sigma, rho);
    }
    
    plhs[0] = mxSigma;
}


void do_update_rho(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mEg(prhs[1]);
    MArray mSigma(prhs[2]);
    MArray mEp(prhs[3]);
    
    matlab_graph_repr egr(mEg);    
    if (egr.weight_class() != mxDOUBLE_CLASS)
    {
        mexErrMsgIdAndTxt("gatrwa:invalidarg", 
                "The edge weights should be of double class.");
    }  
        
    graph_t eg = egr.to_cref_wadjlist_ud<double>();
    
    const double *sigma = mSigma.get_data<double>();
    const double *ep = mEp.get_data<double>();
    
    int m = (int)num_edges(eg);
    
    mxArray *mxRho = mxCreateDoubleMatrix(m, 1, mxREAL);
    double *rho = mxGetPr(mxRho);
    for (int i = 0; i < m; ++i)
    {
        edge_t e(i);
        rho[i] = compute_rho(eg, e, sigma, ep);
    }
    
    plhs[0] = mxRho;
}


void do_comb_update(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mEg(prhs[1]);
    MArray mJdv(prhs[2]);
    mxArray *mxSigma = mxDuplicateArray(prhs[3]);
    mxArray *mxRho = mxDuplicateArray(prhs[4]);
    MArray mEp(prhs[5]);
    MArray mOrd(prhs[6]);
    
    matlab_graph_repr egr(mEg);    
    if (egr.weight_class() != mxDOUBLE_CLASS)
    {
        mexErrMsgIdAndTxt("gatrwa:invalidarg", 
                "The edge weights should be of double class.");
    }          
    graph_t eg = egr.to_cref_wadjlist_ud<double>();
    
    const double *Jdv = mJdv.get_data<double>();
    double *sigma = mxGetPr(mxSigma);
    double *rho = mxGetPr(mxRho);
    const double *ep = mEp.get_data<double>();
    const int *ord = mOrd.get_data<int>();
    
    int N = mOrd.nelems();
    int m = (int)num_edges(eg);
    
    // compute
    
    for (int i = 0; i < N; ++i)
    {
        vertex_t v(i);
        
        sigma[i] = compute_sigma(eg, m, v, Jdv, sigma, rho);
        
        graph_t::out_edge_iterator eep, eend;
        for (boost::tie(eep, eend) = out_edges(v, eg); eep != eend; ++eep)
        {
            edge_t e = *eep;
            int ei = e.i < m ? e.i : e.i - m;
            
            rho[ei] = compute_rho(eg, ei, sigma, ep);
        }        
    }
    
    plhs[0] = mxSigma;
    plhs[1] = mxRho;
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
 *
 * If code == 2: update rho(s)
 *
 *  Input
 *      [1]: eg:        obj.egraph
 *      [2]: sigma:     obj.sigma
 *      [3]: ep:        obj.eprob
 *  Output
 *      [0]: rho:       updated rho vector
 *
 * 
 * If code == 3: combined update of sigma(s) and rho(s)
 *
 * Input
 *     [1]: eg:         obj.egraph
 *     [2]: Jdv:        obj.Jdv
 *     [3]: sigma:      obj.sigma
 *     [4]: rho:        obj.rho
 *     [5]: ep:         obj.eprob
 *     [6]: ord:        order of updating
 * Output
 *     [0]: sigma:      updated sigma vector
 *     [1]: rho:        updated rho vector
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
    else if (code == 2)
    {
        do_update_rho(nlhs, plhs, nrhs, prhs);
    }
    else if (code == 3)
    {
        do_comb_update(nlhs, plhs, nrhs, prhs);
    }
}




