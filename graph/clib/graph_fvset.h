/********************************************************************
 *  
 *  graph_fvset.h
 *
 *  The Feedback vertex set algorithm
 *
 *  Created by Dahua Lin, on Nov 4, 2010
 *
 ********************************************************************/

#ifndef SMI_FEEDBACK_VERTEX_SET_H
#define SMI_FEEDBACK_VERTEX_SET_H

// #define SMI_FVSET_MATLAB_DEBUG


#include "graph_base.h"
#include "../../base/clib/heaps.h"

#include <functional>
#include <vector>
#include <valarray>
#include <queue>
#include <cmath>

namespace smi
{        

    
/**
 * The dynamically pruned graph, which allows one remove a vertex
 * from it at a time, and maintains a heap with max-score.
 */
template<typename TGraph, typename TComp>
class DynamicScorePruneGraph
{
public:
    typedef TGraph graph_type;
    typedef typename TGraph::edge_weight_type score_type;
    typedef graph_size_t size_type;        
    typedef BinaryHeap<score_type, std::greater<score_type> > heap_type;
    
private:         
    const graph_type& m_graph;              // the internal graph representation
    
    size_type m_nr;                 // the number of remaining nodes
    std::valarray<bool> m_ingraph;  // the indicator of whether an vertex is in graph                
    std::valarray<size_type> m_rdegs;       // the remaining degrees        
    
    heap_type m_heap;    
    TComp m_scomp;
    
public:
    DynamicScorePruneGraph(const graph_type& g, TComp scomp) 
    : m_graph(g)
    , m_nr(num_vertices(g)), m_ingraph(true, m_nr), m_rdegs(m_nr)
    , m_heap(m_nr), m_scomp(scomp)
    {
        size_type n = m_nr;
        for (graph_size_t i = 0; i < n; ++i)
        {
            m_rdegs[i] = out_degree(vertex_t(i), g);
        }
    }
    
    // basic information
    
    size_type num_nodes() const
    {
        return m_nr;
    }
    
    bool in_graph(vertex_t v) const
    {
        return m_ingraph[v.i];
    }
    
    size_type degree(vertex_t v) const
    {
        return m_rdegs[v.i];
    }
    
    bool empty() const
    {
        return m_nr == 0;
    }
    
    vertex_t max_score_vertex() const
    {
        return m_heap.root_index();
    }    
    
    score_type max_score() const
    {
        return m_heap.root_key();
    }
    
    
    // graph manipulation
        
    void remove_vertex(vertex_t v)
    {
        if (in_graph(v))
        {
            #ifdef SMI_FVSET_MATLAB_DEBUG
                mexPrintf("remove %d\n", v.i+1);
            #endif            
                                
            // remove the vertex v from heap
            if (v.i == m_heap.root_index())
            {
                m_heap.delete_root();
            }
            else
            {
                m_heap.set_key(v.i, m_heap.root_key() + 1);
                m_heap.delete_root();
            }
            
            // mark it as deleted
            
            m_ingraph[v.i] = false;
            -- m_nr;                        
            
            // reduce degree & update score for all its remaining neighbors
            
            typename graph_type::adjacency_iterator ap, aend;
            for (boost::tie(ap, aend) = adjacent_vertices(v, m_graph); ap != aend; ++ap)
            {
                vertex_t c = *ap;
                if (in_graph(c)) 
                {
                    -- m_rdegs[c.i];
                    update_score(c);
                }
            }
        }
    }
        
    
    // the vertices specified by [vbegin, vend) should have been removed
    // this function going to examine their neighbors recursively
    template<typename TIter>
    void clean_graph_from(TIter vbegin, TIter vend)
    {                        
        std::queue<vertex_t> Q;
        
        for (TIter vp = vbegin; vp != vend; ++vp)
        {
            Q.push(*vp);
        }
        
        while (!Q.empty())
        {
            vertex_t v = Q.front(); 
            Q.pop();  
            
            std::vector<vertex_t> cs;
            cs.reserve(out_degree(v, m_graph));
            
            #ifdef SMI_FVSET_MATLAB_DEBUG
                mexPrintf("pop %d\n", v.i+1);
            #endif
            
            typename graph_type::adjacency_iterator ap, aend;
            for (boost::tie(ap, aend) = adjacent_vertices(v, m_graph); ap != aend; ++ap)
            {
                vertex_t c = *ap;
                if (in_graph(c) && degree(c) < 2)
                {
                    cs.push_back(c);
                    Q.push(c);
                    
                    #ifdef SMI_FVSET_MATLAB_DEBUG
                        mexPrintf("push %d\n", c.i+1);
                    #endif
                }
            }
            
            if (!cs.empty())
            {
                for (std::vector<vertex_t>::const_iterator it = cs.begin(); 
                    it != cs.end(); ++it)
                {
                    remove_vertex(*it);
                }
            }            
        }
    }        
                    
    
    void clean_graph()
    {
        typename graph_type::vertex_iterator vp, vend;
        std::vector<vertex_t> leafs;
        leafs.reserve(num_nodes());
        
        for (boost::tie(vp, vend) = vertices(m_graph); vp != vend; ++vp) 
        {
            vertex_t v = *vp;
            if (degree(v) < 2) {
                leafs.push_back(v);
            }
        }
        
        for (std::vector<vertex_t>::const_iterator it = leafs.begin(); it != leafs.end(); ++it) 
        {
            remove_vertex(*it);
        }
        
        clean_graph_from(leafs.begin(), leafs.end());
    }
    
    
    // score management
    
    score_type compute_score(vertex_t v) const
    {
        return m_scomp(m_graph, &(m_ingraph[0]), v);
    }
    
    void initialize_scores()
    {
        typename graph_type::vertex_iterator vp, vend;
                
        size_type n = num_vertices(m_graph);
        std::valarray<score_type> scores(n);        
        for (boost::tie(vp, vend) = vertices(m_graph); !(vp == vend); ++vp)
        {
            vertex_t v = *vp;
            scores[v.i] = compute_score(v);            
        }
                        
        const score_type *sbegin = &(scores[0]);
        const score_type *send = sbegin + n;        
        m_heap.make_heap(sbegin, send);
    }
                
    void update_score(vertex_t v)
    {
        if (in_graph(v))
        {
            score_type s = compute_score(v);
            m_heap.set_key(v.i, s);
        }        
    }
                
    
}; // end class DynamicScorePruneGraph



// score computers

template<typename TGraph>
struct fvs_deg_computer
{
    typedef TGraph graph_type;
    typedef typename TGraph::edge_weight_type score_type;
    
    score_type operator() (const graph_type& g, const bool *in_graph, vertex_t v) const
    {
        typename graph_type::out_edge_iterator ep, eend;
        score_type s = score_type(0);
        
        for (boost::tie(ep, eend) = out_edges(v, g); ep != eend; ++ep)
        {                        
            edge_t e = *ep;
            vertex_t c = target(e, g);
            
            if (in_graph[c.i])
            {
                s += std::abs(g.get_weight(e));
            }
        }
        
        return s;
    }
};



/**
 * select feedback vertices from a graph
 *
 * @param g the input graph 
 * @param scomp the functor to compute vertex score
 * @param nmax the maximum number of feedback vertices to extract
 * @param out_it the output iterator
 *
 * @return the actual number of feedback vertices that have been extracted
 */
template<typename TGraph, typename TComp, typename TVisitor>
graph_size_t select_feedback_vertex(const TGraph& g, 
        TComp scomp, smi::graph_size_t nmax, TVisitor& vis)
{
    // initialization
    
    #ifdef SMI_FVSET_MATLAB_DEBUG            
        mexPrintf("initialization ......\n");
    #endif
        
    DynamicScorePruneGraph<TGraph, TComp> dspg(g, scomp);    
    dspg.initialize_scores();        
    dspg.clean_graph();
            
    // main loop
    
    #ifdef SMI_FVSET_MATLAB_DEBUG            
        mexPrintf("main loop ......\n");
    #endif
    
    smi::graph_size_t n = 0;

    while (n < nmax && !dspg.empty())
    {                
        vertex_t v = dspg.max_score_vertex();
        vis.select_vertex(v, dspg.max_score());
        ++n;
        
        #ifdef SMI_FVSET_MATLAB_DEBUG
            mexPrintf("select %d (score = %g)\n", v.i+1, dspg.max_score());
        #endif                
        
        dspg.remove_vertex(v);
        dspg.clean_graph_from(&v, (&v)+1);                
    }
    
    return n;    
}


}


#endif




