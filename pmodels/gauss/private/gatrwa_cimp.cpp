/********************************************************************
 *
 *  gatrwa_cimp.cpp
 *
 *  The C++ mex implementation of gatrwa_cimp.m
 *
 *  Created by Dahua Lin, on Feb 6, 2011
 *
 ********************************************************************/


#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/matlab/mgraph.h>
#include <bcslib/graph/bgl_port.h>

#include <cmath>

using namespace bcs;
using namespace bcs::matlab;


typedef gr_wadjlist<double, gr_undirected> graph_t;


inline double calc_b(vertex_t v, const graph_t& eg, gr_index_t m, 
        const double *sigma, const double *rho)
{
    double b = 0;
    graph_t::adj_edge_iterator ep, eend;
    for (rbind(ep, eend) = out_edges(v, eg); ep != eend; ++ep)
    {
        edge_t e = *ep;
        double w = eg.weight_of(e);
        
        gr_index_t ei = e.index < m ? e.index : e.index - m;
        double r = rho[ei];
        
        vertex_t vt = target(e, eg);
        double s = sigma[vt.index];
        
        b += w * r * s;
    }
    
    return b;
}


inline double compute_sigma(const graph_t& eg, gr_index_t m, vertex_t v, 
        const double *Jdv, const double *sigma, const double *rho)
{
    double Jv = Jdv[v.index];    
    double b = calc_b(v, eg, m, sigma, rho);        
    double sig = (std::sqrt(b*b + 4 * Jv) - b) / (2 * Jv);
    
    return sig;
}

inline double compute_rho(const graph_t& eg, edge_t e, 
        const double *sigma, const double *ep)
{
    vertex_t s = source(e, eg);
    vertex_t t = target(e, eg);
    
    double w = eg.weight_of(e);
    double a = w * sigma[s.index] * sigma[t.index];
    
    double beta = ep[e.index];
    
    double r = (beta - std::sqrt(beta*beta + 4 * a*a)) / (2 * a);
    return r;
}




void do_compute_bvec(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_mgraph mEg(prhs[1]);
    const_marray mSigma(prhs[2]);
    const_marray mRho(prhs[3]);
        
    if (mEg.weight_type() != mxDOUBLE_CLASS)
    {
        throw mexception("gatrwa:invalidarg", 
            "The edge weights should be of double class.");
    }  
    
    graph_t eg = to_gr_wadjlist<double, gr_undirected>(mEg);
    gr_size_t n = num_vertices(eg);
    gr_size_t m = num_edges(eg);
    
    const double *sigma = mSigma.data<double>();
    const double *rho = mRho.data<double>();
    
    // compute
    
    marray mB = create_marray<double>(n, 1);
    double *b = mB.data<double>();
    
    for (gr_size_t i = 0; i < n; ++i)
    {
        b[i] = calc_b(i, eg, (gr_index_t)m, sigma, rho);
    }    
    
    plhs[0] = mB.mx_ptr();
    
}


void do_update_sigma(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_mgraph mEg(prhs[1]);
    const_marray mJdv(prhs[2]);
    marray mSigma = duplicate(const_marray(prhs[3]));
    const_marray mRho(prhs[4]);
    const_marray mOrd(prhs[5]);
        
    if (mEg.weight_type() != mxDOUBLE_CLASS)
    {
        throw mexception("gatrwa:invalidarg", 
            "The edge weights should be of double class.");
    }  
    
    
    graph_t eg = to_gr_wadjlist<double, gr_undirected>(mEg);
    
    const double *Jdv = mJdv.data<double>();
    double *sigma = mSigma.data<double>();
    const double *rho = mRho.data<double>();
    const int32_t *ord = mOrd.data<int32_t>();
    
    size_t N = mOrd.nelems();
    gr_size_t m = num_edges(eg);
    
    // compute
    
    for (size_t i = 0; i < N; ++i)
    {
        vertex_t v((gr_index_t)ord[i]);        
        sigma[i] = compute_sigma(eg, (gr_index_t)m, v, Jdv, sigma, rho);
    }
    
    plhs[0] = mSigma.mx_ptr();
}


void do_update_rho(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_mgraph mEg(prhs[1]);
    const_marray mSigma(prhs[2]);
    const_marray mEp(prhs[3]);
       
    if (mEg.weight_type() != mxDOUBLE_CLASS)
    {
        throw mexception("gatrwa:invalidarg", 
            "The edge weights should be of double class.");
    }  
        
    graph_t eg = to_gr_wadjlist<double, gr_undirected>(mEg);
    
    const double *sigma = mSigma.data<double>();
    const double *ep = mEp.data<double>();
    
    gr_size_t m = num_edges(eg);
    
    marray mRho = create_marray<double>(m, 1);
    double *rho = mRho.data<double>();
    for (gr_index_t i = 0; i < (gr_index_t)m; ++i)
    {
        edge_t e(i);
        rho[i] = compute_rho(eg, e, sigma, ep);
    }
    
    plhs[0] = mRho.mx_ptr();
}


void do_comb_update(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_mgraph mEg(prhs[1]);
    const_marray mJdv(prhs[2]);
    marray mSigma = duplicate(const_marray(prhs[3]));
    marray mRho = duplicate(const_marray(prhs[4]));
    const_marray mEp(prhs[5]);
    const_marray mOrd(prhs[6]);
       
    if (mEg.weight_type() != mxDOUBLE_CLASS)
    {
        throw mexception("gatrwa:invalidarg", 
            "The edge weights should be of double class.");
    }          
    graph_t eg = to_gr_wadjlist<double, gr_undirected>(mEg);
    
    const double *Jdv = mJdv.data<double>();
    double *sigma = mSigma.data<double>();
    double *rho = mRho.data<double>();
    const double *ep = mEp.data<double>();
    const int32_t *ord = mOrd.data<int32_t>();
    
    size_t N = mOrd.nelems();
    gr_size_t m = num_edges(eg);
    
    // compute
    
    for (size_t i = 0; i < N; ++i)
    {
        vertex_t v((gr_index_t)i);
        
        sigma[i] = compute_sigma(eg, m, v, Jdv, sigma, rho);
        
        graph_t::adj_edge_iterator eep, eend;
        for (rbind(eep, eend) = out_edges(v, eg); eep != eend; ++eep)
        {
            edge_t e = *eep;
            gr_index_t ei = e.index < m ? e.index : e.index - m;
            
            rho[ei] = compute_rho(eg, ei, sigma, ep);
        }        
    }
    
    plhs[0] = mSigma.mx_ptr();
    plhs[1] = mRho.mx_ptr();
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
    int code = (int)(const_marray(prhs[0]).get_scalar<double>());    
    
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




