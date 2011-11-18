/**********************************************************
 *
 *  mincut_kolmogorov_cimp.cpp
 *
 *  The mex wrapper of Boykov's Max-flow Min-cut algorithm
 *
 *  Created by Dahua Lin, on Sep 26, 2010
 *
 **********************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/matlab/mgraph.h>
#include <bcslib/graph/bgl_port.h>

#include "kolmogorov_maxflow.h"

using namespace bcs;
using namespace bcs::matlab;


/**
 * The minimum-cut solver based on Vladimir Kolmogorov's 
 * implementation (version 3.01)
 * 
 * http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software.html
 */
template<class TWeight>
TWeight mincut_kolmogorov(
        const gr_wadjlist<TWeight, gr_undirected>& g, 
        const TWeight *src_weights, 
        const TWeight *sink_weights,
        bool *results)  // false -> source, true -> sink
{
    gr_size_t n = num_vertices(g);
    gr_size_t m = num_edges(g);               

    // create graph in Kolmogorov's form

    int nmax = (int)(n + 2);
    int mmax = (int)(2 * (n + m + 1));

    typedef vkolmogorov::Graph<TWeight, TWeight, TWeight> VkGraph;

    VkGraph g_vk(nmax, mmax);  

    // add nodes
    g_vk.add_node((int)n);

    // add t-links
    for (gr_size_t i = 0; i < n; ++i)
    {
        g_vk.add_tweights((int)i, src_weights[i], sink_weights[i]);
    }

    // add n-links
    
    for (gr_index_t ei = 0; ei < (gr_index_t)m; ++ei)
    {
        edge_t e(ei);

        vertex_t s = source(e, g);
        vertex_t t = target(e, g);

        int si = (int)s.index;
        int ti = (int)t.index;
        TWeight w = g.weight_of(e);

        g_vk.add_edge(si, ti, w, w);            
    }

    // solve the maxflow problem

    TWeight mfv = g_vk.maxflow();

    // extract the results

    for (gr_size_t i = 0; i < n; ++i)
    {
        results[i] = (g_vk.what_segment((int)i) == VkGraph::SINK);
    }

    return mfv;        
}



template<typename TWeight>
void main_delegate(const_mgraph mG, const_marray mSrcWs, const_marray mSinkWs, 
        int nlhs, mxArray *plhs[])
{
    gr_wadjlist<TWeight, gr_undirected> g = to_gr_wadjlist<TWeight, gr_undirected>(mG);    
    gr_size_t n = num_vertices(g);
    
    const TWeight *src_ws = mSrcWs.data<TWeight>();
    const TWeight *sink_ws = mSinkWs.data<TWeight>();
    
    marray mR = create_marray<bool>(1, n);
    bool *results = mR.data<bool>();
    
    TWeight mfv = mincut_kolmogorov(g, src_ws, sink_ws, results);
    
    plhs[0] = mR.mx_ptr();
    if (nlhs >= 2)
        plhs[1] = create_mscalar<TWeight>(mfv).mx_ptr();
}




/************
 *
 *  main entry:
 *
 *  Input:
 *    [0] g:        the graph (in form of edgelist or adjlist, undirected)
 *    [1] src_ws:   the weights to source
 *    [2] sink_ws:  the weights to sink
 *     
 *  Output
 *    [0] results:  the cut results (1 x n bool vector: 0 - source, 1 - sink)
 *    [1] maxflow:  the value of maximum flow (minimum cut)
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{
    // take input
    
    const_mgraph mG(prhs[0]);        
    const_marray mSrcWs(prhs[1]);
    const_marray mSinkWs(prhs[2]);           
                
    // delegate
    
    switch (mG.weight_type())
    {
        case mxDOUBLE_CLASS:
            main_delegate<double>(mG, mSrcWs, mSinkWs, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            main_delegate<float>(mG, mSrcWs, mSinkWs, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            main_delegate<int>(mG, mSrcWs, mSinkWs, nlhs, plhs);
            break;
            
        case mxUINT32_CLASS:
            main_delegate<unsigned int>(mG, mSrcWs, mSinkWs, nlhs, plhs);
            break;
            
        default:
            throw mexception("mincut_kolmogorov:invalidarg", 
                "The weight value should be double, single, int32, or uint32.");
    }
}

BCSMEX_MAINDEF




