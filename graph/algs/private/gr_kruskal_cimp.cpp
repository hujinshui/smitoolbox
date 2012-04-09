/***************************************************************
 *
 *  gr_kruskal_cimp.cpp
 *
 *  The C++ mex implementation of gr_dijkstra_cimp
 *
 *  Created by Dahua Lin, on Jan 26, 2012
 *
 ***************************************************************/

#include "../../clib/smi_graph_mex.h"
#include <bcslib/graph/graph_minimum_span_trees.h>
#include <bcslib/array/amap.h>
#include <bcslib/base/smart_ptr.h>
#include <vector>

using namespace bcs;
using namespace bcs::matlab;
using namespace smi;


// monitor agents

class krus_monitor
{    
public:
    krus_monitor(std::vector<int32_t>& eds) 
    : edges(eds)
    {
    }
    
    bool examine_edge(const vertex_t&, const vertex_t&, const edge_t&) 
    {
        return true;
    }
    
    bool add_edge(const vertex_t&, const vertex_t&, const edge_t& e)
    {
        edges.push_back(e.id);
        return true;
    }
    
private:
    std::vector<int32_t>& edges;
};

class krus_monitor2
{
public:
    krus_monitor2(std::vector<int32_t>& eds, 
            disjoint_set_forest<vertex_t>& ds, int K_) 
    : edges(eds), dsets(ds), K(K_)
    {
    }
    
    bool examine_edge(const vertex_t&, const vertex_t&, const edge_t&) 
    {
        return true;
    }
    
    bool add_edge(const vertex_t&, const vertex_t&, const edge_t& e)
    {
        edges.push_back(e.id);
        return dsets.ncomponents() > K;
    }
    
private:
    std::vector<int32_t>& edges;
    disjoint_set_forest<vertex_t>& dsets;
    size_t K;
};


marray get_ccs(disjoint_set_forest<vertex_t>& dsets)
{
    typedef std::vector<int32_t> compvec_t;    
    std::vector<bcs::shared_ptr<compvec_t> > comps;
    
    index_t n = (index_t)dsets.size();
    
    // scan clusters
    
    array1d<int32_t> L(n); 
    mem<int32_t>::zero((size_t)n, L.pbase());
    
    int32_t m = 0;
    
    for (index_t i = 0; i < n; ++i)
    {
        vertex_t v;
        v.id = i + 1;
        
        index_t r = dsets.find_root(v);
        int32_t k = L(r);
        
        if (k == 0)
        {
            k = (L(r) = ++m);
            comps.push_back(bcs::shared_ptr<compvec_t>(new compvec_t()));
        }
        
        comps[k-1]->push_back(int32_t(i+1));        
    }
    
    // convert the results to marray
    
    marray mCCS = create_mcell_array(1, m);
    
    for (int32_t k = 0; k < m; ++k)
    {
        const compvec_t& cvec = *(comps[k]);
        mCCS.set_cell(k, to_matlab_row(cvec));
    }
    
    return mCCS;
}



// core function

template<typename T>
void do_kruskal(const gedgelist_t& G, const T* w, int K, int nlhs, mxArray *plhs[])
{    
    // prepare edge-weight map
    
    caview_map<edge_t, T> ewmap(w, 2 * G.nedges());
    
    // prepare outputs
    
    disjoint_set_forest<vertex_t> dsets(G.nvertices());    
    std::vector<int32_t> edges;
    
    // run
    
    if (K == 1)
    {
        krus_monitor mon(edges);
        kruskal_minimum_span_tree_ex(G, ewmap, dsets, mon);
    }
    else
    {
        krus_monitor2 mon(edges, dsets, K);
        kruskal_minimum_span_tree_ex(G, ewmap, dsets, mon);
    }
    
    // extract outputs
    
    plhs[0] = to_matlab_row(edges).mx_ptr();
    if (nlhs > 1)
    {
        plhs[1] = get_ccs(dsets).mx_ptr();
    }
}



/**
 * The main entry
 *
 * Input
 *   [0] G:         the graph struct
 *   [1] w:         the edge weights
 *   [2] K:         the number of components
 *  
 * Output
 *   [0] edges:     the edges of the spanning tree
 *   [1] ccs:       the cell array of connected components
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mG(prhs[0]);
    gedgelist_t G = mx_to_gedgelist(mG);
    
    const_marray mW(prhs[1]);
    const_marray mK(prhs[2]);
    
    int K = (int)mK.get_scalar<int32_t>();
    
    // solve
    
    mxClassID cid = mW.class_id();
    
    if (cid == mxDOUBLE_CLASS)
    {
        do_kruskal(G, mW.data<double>(), K, nlhs, plhs);
    }
    else if (cid == mxSINGLE_CLASS)
    {
        do_kruskal(G, mW.data<float>(), K, nlhs, plhs);
    }
}

BCSMEX_MAINDEF
        
