/********************************************************************
 *
 *  gr_kruskal_cimp.cpp
 *
 *  The C++ mex implementation for Kruskal's minimum spanning tree
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/matlab/mgraph.h>
#include <bcslib/graph/bgl_port.h>

#include <vector>
#include <iterator>

#include <boost/graph/kruskal_min_spanning_tree.hpp>

using namespace bcs;
using namespace bcs::matlab;

typedef boost::default_color_type color_t;


template<typename TWeight>
marray do_kruskal(const_mgraph mG)
{
    gr_wadjlist<TWeight, gr_undirected> g = to_gr_wadjlist<TWeight, gr_undirected>(mG);
        
    std::vector<edge_t> mst_edges;
    mst_edges.reserve(num_vertices(g));
    
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(mst_edges));
    
    size_t m = mst_edges.size();
    marray mEdges = create_marray<int32_t>(1, m);
    int32_t *eds = mEdges.data<int32_t>();
    
    for (size_t i = 0; i < m; ++i)
    {
        eds[i] = (int32_t)mst_edges[i].index + 1;
    }
    
    return mEdges;
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
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    const_mgraph mG(prhs[0]);
    
    marray mEdges;
    
    switch (mG.weight_type())
    {
        case mxDOUBLE_CLASS:
            mEdges = do_kruskal<double>(mG);
            break;
            
        case mxSINGLE_CLASS:
            mEdges = do_kruskal<float>(mG);
            break;
            
        case mxINT32_CLASS:
            mEdges = do_kruskal<int32_t>(mG);
            break;
            
        case mxUINT32_CLASS:
            mEdges = do_kruskal<uint32_t>(mG);
            break;
            
        default:
            throw mexception("gr_kruskal_mst:invalidarg", 
                "The weight value should be double, single, int32, or uint32.");
    }
    
    plhs[0] = mEdges.mx_ptr();
}

BCSMEX_MAINDEF





