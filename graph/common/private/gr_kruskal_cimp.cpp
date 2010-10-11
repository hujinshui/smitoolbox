/********************************************************************
 *
 *  gr_kruskal_cimp.cpp
 *
 *  The C++ mex implementation for Kruskal's minimum spanning tree
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 ********************************************************************/

#include "../../clib/graph_mex.h"

#include <vector>
#include <valarray>
#include <iterator>

#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace smi;

typedef boost::default_color_type color_t;


template<typename TWeight>
void main_delegate(const matlab_graph_repr& gr, int nlhs, mxArray *plhs[])
{
    CRefAdjList<TWeight, boost::undirected_tag> g = gr.to_cref_wadjlist_ud<TWeight>();
        
    graph_size_t n = num_vertices(g);
    
    std::valarray<int> ranks(n);
    std::valarray<vertex_t> preds(n);        
    
    VertexRefMap<int> rmap = &(ranks[0]);
    VertexRefMap<vertex_t> pmap = &(preds[0]);
        
    using boost::predecessor_map;
    using boost::rank_map;
    
    std::vector<edge_t> mst_edges;
    mst_edges.reserve(n);
    
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(mst_edges), 
            rank_map(rmap).predecessor_map(pmap));
    
    plhs[0] = iter_to_matlab_column(mst_edges.begin(), mst_edges.size(),
            unary_chain(edge_to_source(g), vertex_to_mindex()));
    
    plhs[1] = iter_to_matlab_column(mst_edges.begin(), mst_edges.size(),
            unary_chain(edge_to_target(g), vertex_to_mindex()));
    
    if (nlhs >= 3)
    {
        plhs[2] = iter_to_matlab_column(mst_edges.begin(), mst_edges.size(),
                edge_to_weight(g));
    }    
}



/***
 * main entry
 *
 * Inputs:
 *  [0]: G:     the undirected graph with adjlist representation
 *
 * Outputs:
 *  [0]: s:     the sources of edges
 *  [1]: t:     the targets of edges
 *  [2]: w:     the weights of edges
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    matlab_graph_repr gr(prhs[0]);
    
    switch (gr.weight_class())
    {
        case mxDOUBLE_CLASS:
            main_delegate<double>(gr, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            main_delegate<float>(gr, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            main_delegate<int>(gr, nlhs, plhs);
            break;
            
        case mxUINT32_CLASS:
            main_delegate<unsigned int>(gr, nlhs, plhs);
            break;
            
        default:
            mexErrMsgIdAndTxt("gr_kruskal_mst:invalidarg", 
                    "The weight value should be double, single, int32, or uint32.");
    }
}


