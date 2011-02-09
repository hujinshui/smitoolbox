/********************************************************************
 *
 *  ugprune_cimp.cpp
 *
 *  The C++ mex implementation for ugprune
 *
 *  Created by Dahua Lin, on Feb 7, 2011
 *
 ********************************************************************/

#include "../../clib/graph_mex.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <valarray>

using namespace smi;

typedef CRefAdjList<no_edge_weight> graph_t;
typedef graph_t::out_edge_iterator edge_iter_t;

struct wentry
{
    int i;
    double w;
    
    wentry(int i_, double w_) : i(i_), w(w_) { }
    
    bool operator < (const wentry& rhs) const
    {
        return w < rhs.w;
    }
    
    bool operator > (const wentry& rhs) const
    {
        return w > rhs.w;
    }
};


inline void retain(std::valarray<bool>& r, int i, int m)
{    
    if (i < m)
    {
        r[i] = true;
    }
    else
    {
        r[i - m] = true;
    }    
}


void do_select_edges(const graph_t& G0, 
        int K, const double *w, std::vector<int>& sel)
{
    int n = num_vertices(G0);
    int m = num_edges(G0);
    
    std::valarray<bool> r(false, m);
    std::vector<wentry> ss;    
        
    // retain by weights
    
    for (int i = 0; i < n; ++i)
    {
        vertex_t v(i);        
        
        int d = out_degree(v, G0);
        
        if (d <= K)  // retain all
        {
            edge_iter_t ep, eend;
            for (boost::tie(ep, eend) = out_edges(v, G0); ep != eend; ++ep)
            {
                retain(r, ep->i, m);
            }
        }
        else    // select first K
        {
            edge_iter_t ep, eend;
            for (boost::tie(ep, eend) = out_edges(v, G0); ep != eend; ++ep)
            {
                edge_t e = *ep;
                double ew = w[e.i];
                
                ss.push_back(wentry(e.i, ew));                         
            }
                        
            std::nth_element(ss.begin(), ss.begin() + K, ss.end(), 
                        std::greater<wentry>()); 
            
            for (int i = 0; i < K; ++i)
            {
                retain(r, ss[i].i, m);
            }
            
            ss.clear();
        }        
    }
    
    // fill the sel vector
            
    for (int i = 0; i < m; ++i)  
    {
        if (r[i]) sel.push_back(i);
    }
    
}



/***
 * main entry:
 *
 * Inputs:
 *    [0]:  G0      the input graph;
 *    [1]:  K:      the max degree
 *    [2]:  w:      the weights
 *
 * Output
 *    [0]:  sel:    the array of selected edge indices
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    graph_t G0 = matlab_graph_repr(prhs[0]).to_cref_adjlist();
    MArray mK(prhs[1]);
    MArray mW(prhs[2]);
    
    int K = (int)mK.get_double_scalar();
    const double *w = mW.get_data<double>();
    
    std::vector<int> sel;
    sel.reserve(std::max(num_vertices(G0) * K, num_edges(G0)));
    
    do_select_edges(G0, K, w, sel);
    
    plhs[0] = iter_to_matlab_column(sel.begin(), (int)sel.size());    
}
