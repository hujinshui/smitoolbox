/********************************************************************
 *
 *  gr_dijkstra_cimp.cpp
 *
 *  The C++ mex implementation for Dijkstra shortest-path solving
 *  algorithm
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/matlab/mgraph.h>
#include <bcslib/graph/bgl_port.h>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <valarray>

using namespace bcs;
using namespace bcs::matlab;


template<typename TWeight, typename TDir>
void do_dijkstra(const_mgraph mG, const vertex_t& s, int nlhs, mxArray *plhs[])
{
    gr_wadjlist<TWeight, TDir> g = to_gr_wadjlist<TWeight, TDir>(mG);
    gr_size_t nv = num_vertices(g);
    
    if (nlhs <= 1)
    {
        using boost::distance_map;

        std::valarray<TWeight> dists(TWeight(0), nv);

        vertex_ref_map<TWeight> dist_map(&(dists[0]));

        boost::dijkstra_shortest_paths(g, vertex_t(0), distance_map(dist_map));
        
        marray mDists = create_marray<TWeight>(1, nv); 
        copy_elements(&(dists[0]), mDists.data<TWeight>(), nv);
        
        plhs[0] = mDists.mx_ptr();
        
    }
    else
    {
        using boost::distance_map;
        using boost::predecessor_map;

        std::valarray<TWeight> dists(TWeight(0), nv);
        std::valarray<vertex_t> preds(nv);

        vertex_ref_map<TWeight> dist_map(&(dists[0]));
        vertex_ref_map<vertex_t> pred_map(&(preds[0]));

        boost::dijkstra_shortest_paths(g, vertex_t(0), 
                distance_map(dist_map).predecessor_map(pred_map));   
        
        marray mDists = create_marray<TWeight>(1, nv); 
        copy_elements(&(dists[0]), mDists.data<TWeight>(), nv);
        
        marray mPreds = create_marray<int32_t>(1, nv);
        int32_t *ps = mPreds.data<int32_t>();
        for (gr_index_t i = 0; i < (gr_index_t)nv; ++i)
        {
            ps[i] = preds[i].index + 1;
        }
        
        plhs[0] = mDists.mx_ptr();
        plhs[1] = mPreds.mx_ptr();
    }    
}



/***
 * main entry
 *
 * Inputs:
 *  [0]: G:     the graph with adjlist representation
 *  [1]: s:     the source node
 *
 * Outputs:
 *  [0]: parents:   parents in the search tree
 *  [1]: dists:     the distances to the source in search tree 
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    // take input
    
    const_mgraph mG(prhs[0]);
    vertex_t s = const_marray(prhs[1]).get_scalar<vertex_t>();    
    char dtype = mG.dtype();
    
    // main delegate
    
    if (dtype == 'd')
    {
        switch (mG.weight_type())
        {
            case mxDOUBLE_CLASS:
                do_dijkstra<double, gr_directed>(mG, s, nlhs, plhs);
                break;

            case mxSINGLE_CLASS:
                do_dijkstra<float, gr_directed>(mG, s, nlhs, plhs);
                break;

            case mxINT32_CLASS:
                do_dijkstra<int32_t, gr_directed>(mG, s, nlhs, plhs);
                break;

            case mxUINT32_CLASS:
                do_dijkstra<uint32_t, gr_directed>(mG, s, nlhs, plhs);
                break;

            default:
                throw mexception("gr_dijkstra:invalidarg", 
                    "The weight value should be double, single, int32, or uint32.");
        }
    }
    else if (dtype == 'u')
    {
        switch (mG.weight_type())
        {
            case mxDOUBLE_CLASS:
                do_dijkstra<double, gr_undirected>(mG, s, nlhs, plhs);
                break;

            case mxSINGLE_CLASS:
                do_dijkstra<float, gr_undirected>(mG, s, nlhs, plhs);
                break;

            case mxINT32_CLASS:
                do_dijkstra<int32_t, gr_undirected>(mG, s, nlhs, plhs);
                break;

            case mxUINT32_CLASS:
                do_dijkstra<uint32_t, gr_undirected>(mG, s, nlhs, plhs);
                break;

            default:
                throw mexception("gr_dijkstra:invalidarg", 
                    "The weight value should be double, single, int32, or uint32.");
        }
    }
        
}


BCSMEX_MAINDEF

        
