/********************************************************************
 *  
 *  feedback_vertex_set.h
 *
 *  The Feedback vertex set algorithm
 *
 *  Created by Dahua Lin, on Nov 4, 2010
 *
 ********************************************************************/

#ifndef SMI_FEEDBACK_VERTEX_SET_H
#define SMI_FEEDBACK_VERTEX_SET_H


#include "../../../graph/clib/graph_base.h"
#include "../../../base/clib/heaps.h"

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
template<typename TGraph>
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
    
public:
    DynamicScorePruneGraph(const graph_type& g) 
    : m_graph(g)
    , m_nr(num_vertices(g)), m_ingraph(true, m_nr), m_rdegs(m_nr)
    , m_heap(m_nr)
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
    
    // graph manipulation
        
    void remove_vertex(vertex_t v)
    {
        if (in_graph(v))
        {
            m_heap.delete_root();
            
            typename graph_type::adjacency_iterator ap, aend;
            for (boost::tie(ap, aend) = adjacent_vertices(v, m_graph); ap != aend; ++ap)
            {
                vertex_t c = *ap;
                if (in_graph(c)) 
                    -- m_rdegs[c.i];
            }
            
            m_ingraph[v.i] = false;
            -- m_nr;
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
            vertex_t v = Q.front(); // v should have been removed
            Q.pop();  
            
            typename graph_type::adjacency_iterator ap, aend;
            for (boost::tie(ap, aend) = adjacent_vertices(v, m_graph); ap != aend; ++ap)
            {
                vertex_t c = *ap;
                if (in_graph(c) && degree(c) < 2)
                {
                    remove_vertex(c);
                    Q.push(c);
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
    
    template<typename TComp>
    score_type compute_score(vertex_t v, TComp scomp) const
    {
        return scomp(m_graph, &(m_ingraph[0]), v);
    }
    
    template<typename TComp>
    void initialize_scores(TComp scomp)
    {
        typename graph_type::vertex_iterator vp, vend;
        
        std::valarray<score_type> scores;        
        for (boost::tie(vp, vend) = vertices(m_graph); !(vp == vend); ++vp)
        {
            vertex_t v = *vp;
            scores[v.i] = compute_score<TComp>(v, scomp);            
        }
        
        const score_type *sbegin = &(scores[0]);
        const score_type *send = sbegin + num_vertices(m_graph);        
        m_heap.make_heap(sbegin, send);
    }
    
        
    template<typename TComp>
    void update_score(vertex_t v, TComp scomp)
    {
        if (in_graph(v))
        {
            score_type s = compute_score<TComp>(v, scomp);
            m_heap.set_key(v.i, s);
        }        
    }
    
    
    template<typename TComp>
    void update_neighbor_scores(vertex_t u, TComp scomp)
    {
        typename graph_type::adjacency_iterator ap, aend;
        for (boost::tie(ap, aend) = adjacent_vertices(u, m_graph); ap != aend; ++ap)
        {
            update_score<TComp>(*ap, scomp);
        }
    }
                
    
}; // end class DynamicScorePruneGraph



// score computers

template<typename TGraph>
struct fvs_deg_computer
{
    typedef TGraph graph_type;
    typedef typename TGraph::edge_weight_type score_type;
    
    score_type operator() (const graph_type& g, const bool *in_graph, vertex_t v)
    {
        typename graph_type::out_edge_iterator ep, eend;
        score_type s = score_type(0);
        
        for (boost::tie(ep, eend) = out_edges(v, g); ep != eend; ++ep)
        {                        
            edge_t e = *ep;
            if (in_graph[v.i])
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
template<typename TGraph, typename TComp, typename TInsertIter>
graph_size_t select_feedback_vertex(const TGraph& g, 
        TComp scomp, smi::graph_size_t nmax, TInsertIter out_it)
{
    // initialization
    
    DynamicScorePruneGraph<TGraph> dspg(g);    
    dspg.initialize_scores(scomp);        
    dspg.clean_graph();
            
    // main loop
    
    smi::graph_size_t n;
    while (n < nmax && !dspg.empty())
    {
        vertex_t v = dspg.max_score_vertex();
        *(out_it++) = v;
        
        dspg.remove_vertex(v);
        dspg.clean_graph_from(&v, (&v)+1);
    }
    
    return n;    
}


}


#endif




