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

#include "graph_adjlist.h"

namespace smi
{
    
    
    
    
    
    /**
     * The minimum-cut solver based on Vladimir Kolmogorov's 
     * implementation (version 3.01)
     * 
     * http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software.html
     */
    template<class TWeight>
    TWeight mincut_kolmogorov(
            const CRefEdgeList<TWeight, boost::undirected_tag>& g, 
            const TWeight *src_weights, 
            const TWeight *sink_weights)
    {
        graph_size_t n = num_vertices(g);
        graph_size_t m = num_edges(g);
        
        
        
        for (edge_iterator_t it = g.edges_begin(); it != g.edges_end(); ++it)
        {
            edge_t e = *it;
            
            vertex_t s = source(e, g);
            vertex_t t = target(e, g);
            
            
        }
    }




}

#endif
