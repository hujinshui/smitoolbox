/********************************************************************
 *
 *  graph_mincut_kolmogorov.h
 *
 *  The functions for graph minimum-cut using Kolmogorov's
 *  implementation.
 *  
 *  Created by Dahua Lin, on Oct 14, 2010
 *
 ********************************************************************/

#ifndef SMI_GRAPH_MINCUT_KOLMOGOROV_H
#define SMI_GRAPH_MINCUT_KOLMOGOROV_H

#include "graph_adjlist.h"
#include "kolmogorov_maxflow.h"


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
            const TWeight *sink_weights,
            bool *results)  // false -> source, true -> sink
    {
        graph_size_t n = num_vertices(g);
        graph_size_t m = num_edges(g);               
        
        // create graph in Kolmogorov's form
        
        int nmax = (int)(n + 2);
        int mmax = (int)(2 * (n + m + 1));
        
        typedef vkolmogorov::Graph<TWeight, TWeight, TWeight> VkGraph;
        
        VkGraph g_vk(nmax, mmax);  
        
        // add nodes
        g_vk.add_node((int)n);
        
        // add t-links
        for (graph_size_t i = 0; i < n; ++i)
        {
            g_vk.add_tweights((int)i, src_weights[i], sink_weights[i]);
        }
        
        // add n-links
        for (edge_iterator_t it = g.edges_begin(); it != g.edges_end(); ++it)
        {
            edge_t e = *it;
            
            vertex_t s = source(e, g);
            vertex_t t = target(e, g);
            
            int si = (int)s.i;
            int ti = (int)t.i;
            TWeight w = g.get_weight(e);
            
            if (si < ti)
            {
                g_vk.add_edge(si, ti, w, w);
            }            
        }
        
        // solve the maxflow problem
        
        TWeight mfv = g_vk.maxflow();
        
        // extract the results
        
        for (graph_size_t i = 0; i < n; ++i)
        {
            results[i] = (g_vk.what_segment((int)i) == VkGraph::SINK);
        }
        
        return mfv;        
    }

}

#endif

