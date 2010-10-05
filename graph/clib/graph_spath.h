/********************************************************************
 *  
 *  graph_spath.h
 *
 *  The header files for graph shortest path algorithms
 *
 *  Created by Dahua Lin, on Oct 4, 2010
 *
 ********************************************************************/

#ifndef SMI_GRAPH_SPATH_H
#define SMI_GRAPH_SPATH_H

#include "graph_search.h"

namespace smi
{

/**
 * The base class for all single source shortest path iterator
 */    
template<typename TWeight>    
class SingleSource_SPathIterator 
{
public:
    SingleSource_SPathIterator(int n) 
    : m_n(n), m_s(-1)
    , m_discovered(n), m_closed(n), m_close_seq(n), m_d(n), m_preds(n)
    {
        m_discovered.set_zeros();
        m_closed.set_zeros();        
    }
    
    int num_nodes() const
    {
        return m_n;
    }
    
    int source_node() const
    {
        return m_s;
    }
    
    int is_discovered(int v) const
    {
        return m_discovered[v];
    }
    
    int is_closed(int v) const
    {
        return m_closed[v];
    }
    
    TWeight distance_of(int v) const
    {
        return m_d[v];
    }    
        
    int predecessor_of(int v) const
    {
        return m_preds[v];
    }
    
    const SeqList<int>& close_sequence() const
    {
        return m_close_seq;
    }
    
    int num_closed() const
    {
        return m_close_seq.size();
    }
    
protected:
    
    void set_source(int s)
    {
        m_s = s;
    }
    
    bool relax(int u, int v, TWeight w)
    {
        TWeight new_d = m_d[u] + w;
        
        if (!is_discovered(v))
        {
            discover_node(u, v, new_d);            
            return true;
        }
        else
        {
            if (new_d < m_d[v])
            {                
                m_d[v] = new_d;
                m_preds[v] = u;
                
                return true;
            }            
            else
            {
                return false;
            }
        }
    }
    
    void discover_node(int pred, int v, TWeight d0)
    {
        m_discovered[v] = true;
        m_preds[v] = pred;
        m_d[v] = d0;
    }    
        
    void close_node(int v)
    {
        m_closed[v] = true;
        m_close_seq.add(v);
    }
                
private:    
    // disable copying
    SingleSource_SPathIterator(const SingleSource_SPathIterator<TWeight>& );
    SingleSource_SPathIterator& operator = (const SingleSource_SPathIterator<TWeight>& );
    
private:
    int m_n;    // the number of nodes
    int m_s;    // the source node index
    
    Array<bool> m_discovered;   // whether a node is discovered (initial bound assigned)
    Array<bool> m_closed;   // whether a node is closed (s-path determined)
    
    SeqList<int> m_close_seq;   // the sequence of closed nodes    
    
    Array<TWeight> m_d;     // the distance bounds
                            // the value is meaningful only after discovered
                            // the value equals shortest distance after closed
    
    Array<int> m_preds;     // the predecessors
        
            
}; // end class SingleSource_SPathIterator
    
    



/**
 * The solving iterator for single-source shortest path of DAG
 */    
template<typename TWeight>
class DAG_SPathIterator : public SingleSource_SPathIterator<TWeight>
{
public:
    explicit DAG_SPathIterator(const WAdjList<TWeight>& adjlist)
    : SingleSource_SPathIterator(adjlist.nnodes())
    , m_adjlist(adjlist), m_tord(adjlist.nnodes())
    {
    }       
    
    const SeqList<int>& topological_order()
    {
        return m_tord;
    }     
    
    bool done() const
    {
        return this->num_closed() == m_tord.size();
    }
    
    
    /**
     * Initialize the object for search from source s
     *
     * @return true if the accessible nodes from s form a DAG,
     *         false otherwise
     *
     * @remark the instance can be re-initialized for a different source.
     */
    bool initialize(int s)
    {
        set_source(s);
        
        DFSIteratorEx dfs_it(m_adjlist);
        dfs_it.set_seed(s);
        
        if (dfs_it.find_loop() >= 0)
            return false;
        
        const SeqList<int>& ford = dfs_it.finish_order();
        int nf = ford.size();
        for (int i = 0; i < nf; ++i)
        {
            m_tord[i] = ford[nf-i-1];
        }
        
        this->discover_node(-1, s, 0); 
        
        return true;
    }
    
    
    /**
     * Move the iterator forward (run next step)
     *
     * @return the node that is closed in this step, or -1 if already done     
     */
    int next()
    {        
        int ti = this->num_closed();
        
        if (ti < m_tord.size())
        {
            int u = m_tord[ti];
        
            int nnb = m_adjlist.neighbor_num(u);
            const int *nbs = m_adjlist.neighbor_nodes(u);
            const TWeight *ws = m_adjlist.neighbor_weights(u);
        
            for (int j = 0; j < nnb; ++j)
            {
                int v = nbs[j];
                TWeight w = ws[j];
                        
                this->relax(u, v, w);
            }
            
            this->close_node(u);
            return u;
        }
        else
        {
            return -1;
        }        
    }
    
    
    /**
     * Solve the shortest path to all accessible nodes     
     */
    void solve_all()
    {
        while (next() >= 0);
    }
    
private:
    DAG_SPathIterator(const DAG_SPathIterator<TWeight>& );
    DAG_SPathIterator& operator = (const DAG_SPathIterator<TWeight>& );
    
private:
    const WAdjList<TWeight>& m_adjlist;      // the adjacency list of the graph
    SeqList<int> m_tord;            // the topological order of the nodes                   
    
}; // end class DAG_SPathIterator
    





};

#endif

