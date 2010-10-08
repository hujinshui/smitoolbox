/********************************************************************
 *
 *  graph_mst.h
 *
 *  The classes for Minimum Spanning Tree
 *
 *  Created by Dahua Lin, on Oct 7, 2010
 *
 ********************************************************************/


#ifndef SMI_GRAPH_MST_H
#define SMI_GRAPH_MST_H


#include "graphs.h"
#include "../../base/clib/array.h"
#include "../../base/clib/heaps.h"


namespace smi
{


template<typename TWeight>
class Prim_MSTIterator
{
public:
    struct entry
    {
        entry() { }
        
        entry(int u, int v, TWeight w) : parent(u), node(v), weight(w) { }
        
        int parent;
        int node;        
        TWeight weight;   
        
        bool operator < (const entry& rhs) const
        {
            return weight < rhs.weight;
        }
    };
    
    typedef WEdge<TWeight> wedge_t;
    typedef BinaryHeap<entry, std::less<entry> > heap_t;
    typedef typename heap_t::index_type index_t;
    
public:
    Prim_MSTIterator(const WAdjList<TWeight>& G)
    : m_adjlist(G), m_root(-1), m_heap(G.nnodes())
    , m_discovered(G.nnodes()), m_included(G.nnodes()), m_imap(G.nnodes())
    {
        m_discovered.set_zeros();
        m_included.set_zeros();
    }
    
    
    int root_node() const
    {
        return m_root;
    }
    
    bool is_discovered(int v) const
    {
        return m_discovered[v];
    }
    
    bool is_included(int v) const
    {
        return m_included[v];
    }
    
    
    /**
     * Initialize the algorithm iterator with an appointed root node
     *
     * Can initialize with new seed when the running on a connected
     * component has been done.
     */
    void initialize(int root)
    {
        if (!is_included(root))
        {
            m_root = root;        
            m_discovered[root] = true;
        
            include_node(root);
        }                
    }
    
    
    /**
     * Locate the next edge in MST.
     *
     * @return true is next edge is found, false otherwise (at the end)
     */
    bool next(wedge_t& edge)
    {
        if (!m_heap.empty())
        {
            const entry& e = m_heap.root_key();
        
            edge.s = e.parent;
            edge.t = e.node;
            edge.w = e.weight;
            
            m_heap.delete_root();
            
            include_node(edge.t);                        
            
            return true;
        }
        else
        {
            return false;
        }        
    }
    
    
    template<typename OutputIter>
    int solve_all(OutputIter it)
    {
        wedge_t edge;
       
        int c = 0;
        while (next(edge))
        {
            *(it++) = edge;
            ++c;
        }
        
        return c;
    }
    
    
private:
    
    void include_node(int u)
    {                
        m_included[u] = true;
        
        int nnb = m_adjlist.neighbor_num(u);
        const int *nbs = m_adjlist.neighbor_nodes(u);
        const TWeight *ws = m_adjlist.neighbor_weights(u);
        
        for (int i = 0; i < nnb; ++i)
        {
            int v = nbs[i];
            if (!is_included(v))
            {
                TWeight w = ws[i];
                
                if (is_discovered(v))
                {
                    index_t idx = m_imap[v];
                    entry& e = m_heap.get_key(idx);
                    if (w < e.weight)
                    {
                        e.parent = u;
                        e.weight = w;
                    }
                    m_heap.notify_update(idx);
                }
                else
                {
                    m_imap[v] = m_heap.add_key(entry(u, v, w));
                    m_discovered[v] = true;
                }
            }
        }
    }
    
    
private:
    const WAdjList<TWeight>& m_adjlist;
    int m_root;
    
    heap_t m_heap;
    Array<bool> m_discovered;
    Array<bool> m_included;
    Array<index_t> m_imap;
    
}; // end class Prim_MSTIterator



};

#endif

