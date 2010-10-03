/********************************************************************
 *  
 *  graph_search.h
 *
 *  The header files for graph search (BFS/DFS) algorithms
 *
 *  Created by Dahua Lin, on Oct 2, 2010
 *
 ********************************************************************/

#ifndef SMI_GRAPH_SEARCH_H
#define SMI_GRAPH_SEARCH_H

#include "graphs.h"
#include "../../../base/clib/array.h"

#define _SECURE_SCL 0
#include <vector>
#include <deque>

namespace smi
{
  
/**
 * The iterator of graph nodes using BFS order
 *
 * @remarks It is a forward-input iterator
 */     
class BFSIterator
{
public:
        
    BFSIterator(const AdjList& adjlist) 
    : m_adjlist(adjlist), m_Q(adjlist.nnodes()), m_closed(adjlist.nnodes())
    {
        m_closed.set_zeros();
    }   
    
    
    /**
     * add a new node as seed, and tag v as closed
     *
     * @param v the index of the node to be added as seed
     * @return whether v is added (true if v is open before adding) 
     *
     * @remarks      
     *  - this function adds v only when v remains open       
     */
    bool add_seed(int v)
    {
        if (!is_closed(v))
        {
            enqueue(v);
        }        
    }
    
    
    /**
     * Get the next node and move the iterator forward
     *
     * @return the index of next node, or -1 when no next node
     */
    int next()
    {
        if (!m_Q.empty())
        {
            int s = m_Q.front();
            m_Q.pop_front();
            
            int c = m_adjlist.neighbor_num(s);
            if (c > 0)
            {
                const int *ts = m_adjlist.neighbor_nodes(s);
                for (int i = 0; i < c; ++i)
                {
                    int t = ts[i];
                    if (!is_closed(t))
                    {
                        enqueue(t);
                    }   
                }
            }
            
            return s;
        }
        else
        {
            return -1;
        }
    }
    
    
public:
    
    bool is_closed(int v) const
    {
        return m_closed[v];
    }
    
    bool is_queue_empty() const
    {
        return m_Q.empty();
    }
        
    int queue_length() const
    {
        return (int)m_Q.size();
    }
    
    
private:
    
    void enqueue(int v)
    {
        m_Q.push_back(v);
        m_closed[v] = true;        
        
        mexPrintf("Q: ");
        for (std::deque<int>::const_iterator it = m_Q.begin(); it != m_Q.end(); ++it)
        {
            mexPrintf("%d ", *it);            
        }
        mexPrintf("\n");
    }
    
                
private:
    const AdjList& m_adjlist;
    std::deque<int> m_Q;
    
    Array<bool> m_closed;
    
}; // end class BFSIterator
    


}


#endif

