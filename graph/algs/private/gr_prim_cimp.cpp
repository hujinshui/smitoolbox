/***************************************************************
 *
 *  gr_prim_cimp.cpp
 *
 *  The C++ mex implementation of gr_prim_cimp
 *
 *  Created by Dahua Lin, on Jan 26, 2012
 *
 ***************************************************************/

#include "../../clib/smi_graph_mex.h"
#include <bcslib/graph/graph_minimum_span_trees.h>
#include <bcslib/array/amap.h>
#include <vector>

using namespace bcs;
using namespace bcs::matlab;
using namespace smi; 


// agents

template<typename T>
class prim_monitor
{
public:
    prim_monitor(std::vector<int32_t>& eds, 
            const bool *vmsk,
            const T* eweights,
            const gint& max_size, 
            const T& max_w)            
    : edges(eds), m_vmsk(vmsk), m_eweights(eweights)
    , m_max_size(max_size), m_max_w(max_w)
    {
    }
    
    bool examine_edge(const vertex_t&, const vertex_t& v, const edge_t& e) 
    {
        return (!m_vmsk || m_vmsk[v.index()]) && 
               (m_max_w < 0 || m_eweights[e.index()] <= m_max_w);
        return true;
    }
    
    bool add_edge(const vertex_t&, const vertex_t&, const edge_t& e)
    {
         
        edges.push_back(e.id);
        return m_max_size < 0 || ((gint)edges.size() + 1 < m_max_size);
    }
    
    std::vector<int32_t>& edges;
    
private:
    const bool *m_vmsk;
    const T* m_eweights;
    gint m_max_size;
    T m_max_w;
};


// core algorithm

template<typename T>
void do_prim(const ginclist_t& G, const T *w, const vertex_t& rv, 
        const bool* vmsk, gint max_size, T max_w, 
        std::vector<int32_t>& edges)
{
    size_t ma = max_size < 0 ? G.nvertices() : max_size;
    edges.reserve(ma);
    
    caview_map<edge_t, T> ewmap(w, G.nedges() * 2);
    
    prim_monitor<T> mon(edges, vmsk, w, max_size, max_w);
    
    prim_minimum_span_tree_ex(G, ewmap, rv, mon);
}



/**
 * The main entry
 *
 * Input
 *   [0] G:         the graph struct
 *   [1] w:         the edge weights
 *   [2] rv:        the root vertex
 *   [3] msk:       the vertex mask
 *   [4] max_size:  the maximum component size
 *   [5] max_w:     the maximum weight
 *  
 * Output
 *   [0] edges:     the edges of the spanning tree
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mG(prhs[0]);
    const_marray mW(prhs[1]);
    const_marray mRv(prhs[2]);
    const_marray mMsk(prhs[3]);
    const_marray mMaxSize(prhs[4]);
    const_marray mMaxW(prhs[5]);
    
    ginclist_t G = mx_to_ginclist(mG);    
    gint rv_id = mRv.get_scalar<gint>();    
    vertex_t rv = make_gvertex(rv_id);
    const bool *vmsk = (mMsk.is_empty() ? (const bool*)(0) : mMsk.data<bool>());    
    
    gint max_size = (gint)mMaxSize.get_scalar<gint>(); 
      
    // solve
    
    mxClassID cid = mW.class_id();
    
    std::vector<int32_t> edges;
    
    if (cid == mxDOUBLE_CLASS)
    {        
        double max_w = mMaxW.get_scalar<double>();
        do_prim(G, mW.data<double>(), rv, vmsk, max_size, max_w, edges);
    }
    else if (cid == mxSINGLE_CLASS)
    {
        float max_w = mMaxW.get_scalar<float>();
        do_prim(G, mW.data<float>(), rv, vmsk, max_size, max_w, edges);
    }
    
    // extract output
    
    plhs[0] = to_matlab_row(edges).mx_ptr();
    
}

BCSMEX_MAINDEF

