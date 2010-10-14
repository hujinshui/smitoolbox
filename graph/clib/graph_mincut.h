/********************************************************************
 *
 *  graph_mincut.h
 *
 *  The functions for graph minimum-cut
 *  
 *  Created by Dahua Lin, on Oct 14, 2010
 *
 ********************************************************************/

#ifndef SMI_GRAPH_MINCUT_H
#define SMI_GRAPH_MINCUT_H

#include "graph_base.h"

namespace smi
{
    
    void mincut_kolmogorov_double()

    /**
     * The minimum-cut solver based on Vladimir Kolmogorov's 
     * implementation (version 3.01)
     * 
     * http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software.html
     */
    template<class TGraph, class VWeightMap>
    TWeight mincut_kolmogorov(
            const TGraph& g, 
            const VWeightMap& source_weights,
            const VWeightMap& sink_weights)
    {
        typedef typename boost::graph_traits<TGraph>::vertices_size_type vsize_t;
        typedef typename boost::graph_traits<TGraph>::edges_size_type esize_t;
        
        typedef typename boost::graph_traits<TGraph>::vertex_descriptor vertex_ty;
        typedef typename boost::graph_traits<TGraph>::edge_descriptor edge_ty;        
        typedef typename boost::graph_traits<TGraph>::vertex_iterator vertex_it;
        typedef typename boost::graph_traits<TGraph>::edge_iterator edge_it;
        
        typedef typename boost::property_map<TGraph, boost::vertex_index_t>::const_type vindexmap_ty;
        vindexmap_ty vindexmap = get(g, boost::vertex_index);
        
        typedef typename boost::property_map<TGraph, boost::edge_weight_t>::const_type eweightmap_ty;
        eweightmap_ty eweightmap = get(g, boost::edge_weight);
        
        vertex_it vit, vend;
        for (boost::tie(vit, vend) = g.vertices; vit != vend; ++ vit)
        {
            
        }
        
        return TWeight(0);
    }




}

#endif
