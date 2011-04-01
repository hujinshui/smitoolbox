/********************************************************************
 *
 *  gr_prim_cimp.cpp
 *
 *  The C++ mex implementation for Prim's minimum spanning tree
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/matlab/mgraph.h>
#include <bcslib/graph/bgl_port.h>

#include <vector>
#include <valarray>
#include <iterator>

#include <boost/graph/prim_minimum_spanning_tree.hpp>

using namespace bcs;
using namespace bcs::matlab;


class PrimVisitor : public boost::default_dijkstra_visitor
{
public:
    
    PrimVisitor(edge_t *eds): m_edges(eds)
    {        
    }
    
    template<class Graph>
    void edge_relaxed(const edge_t& e, const Graph& g)
    {        
        vertex_t t = target(e, g);
        m_edges[t.index] = e;
    }
        
private:
    edge_t *m_edges;
};


template<typename TWeight>
void do_prim(const_mgraph mG, const vertex_t& root, int nlhs, mxArray *plhs[])
{
    gr_wadjlist<TWeight, gr_undirected> g = to_gr_wadjlist<TWeight, gr_undirected>(mG);
    gr_size_t nv = num_vertices(g);
    
    std::valarray<vertex_t> preds(nv);
    
    if (nlhs <= 1)
    {
        boost::prim_minimum_spanning_tree(g, vertex_ref_map<vertex_t>(&(preds[0])) );
    }
    else
    {
        std::valarray<edge_t> erec(nv);      
        gr_index_t ne = (gr_index_t)g.nedges();
        
        boost::prim_minimum_spanning_tree(g, vertex_ref_map<vertex_t>(&(preds[0])),
                boost::visitor(PrimVisitor(&(erec[0]))) );
        
        marray mEdgeMap = create_marray<int32_t>(1, nv);
        int32_t *eds = mEdgeMap.data<int32_t>();
        for (gr_size_t i = 0; i < nv; ++i)
        {
            gr_index_t ei = erec[i].index + 1;
            if (ei >= ne)
            {
                ei -= ne;
            }
            eds[i] = ei; 
        }
        
        plhs[1] = mEdgeMap.mx_ptr();
    }
    
    marray mPredMap = create_marray<int32_t>(1, nv);
    int32_t *ps = mPredMap.data<int32_t>();
    for (gr_size_t i = 0; i < nv; ++i)
    {
        ps[i] = preds[i].index + 1;
    }
    
    plhs[0] = mPredMap.mx_ptr();
}



/***
 * main entry
 *
 * Inputs:
 *  [0]: G:     the undirected graph with adjlist representation
 *  [1]: r:     the root node of the MST
 *
 * Outputs:
 *  [0]: preds:  the map of preceding vertices
 *  [1]: eds:    the map of preceding edges
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    const_mgraph mG(prhs[0]);
    vertex_t root = (gr_index_t)(const_marray(prhs[1]).get_scalar<int32_t>());
    
    switch (mG.weight_type())
    {
        case mxDOUBLE_CLASS:
            do_prim<double>(mG, root, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            do_prim<float>(mG, root, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            do_prim<int32_t>(mG, root, nlhs, plhs);
            break;
            
        case mxUINT32_CLASS:
            do_prim<uint32_t>(mG, root, nlhs, plhs);
            break;
            
        default:
            throw mexception("gr_prim_mst:invalidarg", 
                "The weight value should be double, single, int32, or uint32.");
    }
}

BCSMEX_MAINDEF




