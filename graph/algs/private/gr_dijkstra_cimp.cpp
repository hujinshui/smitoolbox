/***************************************************************
 *
 *  gr_dijkstra_cimp.cpp
 *
 *  The C++ mex implementation of gr_dijkstra_cimp
 *
 *  Created by Dahua Lin, on Jan 26, 2012
 *
 ***************************************************************/

#include "../../clib/smi_graph_mex.h"
#include <bcslib/graph/graph_shortest_paths.h>
#include <bcslib/array/amap.h>
#include <limits>
#include <vector>

using namespace bcs;
using namespace bcs::matlab;
using namespace smi;

template<typename T>
struct dijks_agentL1 : public trivial_dijkstra_agent<ginclist_t, T>
{
};

template<typename T>
struct dijks_agentL2 : public trivial_dijkstra_agent<ginclist_t, T>
{    
    std::vector<gint>& vs;
    
    dijks_agentL2(gint n, std::vector<gint>& vs_) : vs(vs_)
    {
        vs.reserve(n);
    }
    
    bool enroll(const vertex_t& u, const T& )
    {
        vs.push_back(u.id);
        return true;
    }
};


template<typename T>
struct dijks_agentL3 : public trivial_dijkstra_agent<ginclist_t, T>
{
    std::vector<gint>& vs;
    std::vector<gint>& preds;
    
    array_map<vertex_t, vertex_t> pred_map;
            
    dijks_agentL3(gint n, std::vector<gint>& vs_, std::vector<gint>& preds_) 
    : vs(vs_), preds(preds_), pred_map(n, make_gvertex(0))
    {
        vs.reserve(n);
    }
    
    bool enroll(const vertex_t& u, const T& )
    {
        vs.push_back(u.id);
        preds.push_back(pred_map[u].id);
        return true;
    }
    
    bool discover(const vertex_t& u, const vertex_t& v, const edge_t&, const T& )
    {
        pred_map[v] = u;
        return true;
    }
    
    bool relax(const vertex_t& u, const vertex_t& v, const edge_t&, const T& )
    {
        pred_map[v] = u;
        return true;
    }
};




template<typename T>
void do_dijkstra(const ginclist_t& g, const T *w, 
        index_t ns, const vertex_t *s, int nlhs, mxArray *plhs[])
{
    gint n = g.nvertices();
    gint m = g.nedges();
    
    // prepare edge-weight map
    
    gint ma = g.is_directed() ? m : 2 * m;
    caview_map<edge_t, T> ewmap(w, ma);
    
    // prepare outputs
    
    marray mLens = create_marray<T>(1, n);
    aview_map<vertex_t, T> plmap(mLens.data<T>(), n);
    
    std::vector<gint> vs; 
    std::vector<gint> preds;
    
    // run Dijkstra's algorithm
    
    T inf = std::numeric_limits<T>::infinity();
    
    if (nlhs <= 1)
    {
        dijks_agentL1<T> agent;
        dijkstra_shortest_paths(g, ewmap, plmap, inf, agent, s, s + ns);
    }
    else if (nlhs == 2)
    {                
        dijks_agentL2<T> agent(n, vs);
        dijkstra_shortest_paths(g, ewmap, plmap, inf, agent, s, s + ns);
    }
    else if (nlhs == 3)
    {
        dijks_agentL3<T> agent(n, vs, preds);
        dijkstra_shortest_paths(g, ewmap, plmap, inf, agent, s, s + ns);
    }
    
    // extract outputs
    
    plhs[0] = mLens.mx_ptr();    
    if (nlhs > 1) plhs[1] = to_matlab_row(vs).mx_ptr();    
    if (nlhs > 2) plhs[2] = to_matlab_row(preds).mx_ptr();

}


/**
 * The main entry
 *
 * Input
 *   [0] G:         the graph struct
 *   [1] w:         the edge weights
 *   [2] s:         tht vectex of sources
 *  
 * Output
 *   [0] lens:      the shortest path lengths to all vertices
 *   [1] vs:        the vector of vertices in ascending order of shortest
 *                  path lengths
 *   [2] preds:     the vector of corresponding predecessors
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mG(prhs[0]);
    ginclist_t G = mx_to_ginclist(mG);
    
    const_marray mW(prhs[1]);
    const_marray mS(prhs[2]);
    
    index_t ns = mS.nelems();
    const vertex_t* s = mS.data<vertex_t>();
    
    // solve
    
    mxClassID cid = mW.class_id();
    
    if (cid == mxDOUBLE_CLASS)
    {
        do_dijkstra(G, mW.data<double>(), ns, s, nlhs, plhs);
    }
    else if (cid == mxSINGLE_CLASS)
    {
        do_dijkstra(G, mW.data<float>(), ns, s, nlhs, plhs);
    }
}

BCSMEX_MAINDEF

        