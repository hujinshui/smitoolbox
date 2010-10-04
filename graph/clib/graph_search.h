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
#include "../../base/clib/array.h"
#include "../../base/clib/data_structs.h"


namespace smi
{
  
/**
 * The iterator of graph nodes using BFS order
 */     
class BFSIterator
{
public:
        
    BFSIterator(const AdjList& adjlist) 
    : m_adjlist(adjlist), m_Q(adjlist.nnodes()), m_discovered(adjlist.nnodes())
    {
        m_discovered.set_zeros();
    }   
    
    
    /**
     * add a new node as seed, and tag v as discovered
     *
     * @param v the index of the node to be added as seed
     * @return whether v is added (true if v is undiscovered before adding) 
     *
     * @remarks      
     *  - this function adds v only when v remains undiscovered      
     */
    bool add_seed(int v)
    {
        if (!is_discovered(v))
        {
            discover_node(v);
            return true;
        } 
        else
        {
            return false;
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
            m_Q.dequeue();
            
            int c = m_adjlist.neighbor_num(s);
            if (c > 0)
            {
                const int *ts = m_adjlist.neighbor_nodes(s);
                for (int i = 0; i < c; ++i)
                {
                    int t = ts[i];
                    if (!is_discovered(t))
                    {
                        discover_node(t);
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
    
    bool is_discovered(int v) const
    {
        return m_discovered[v];
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
    
    void discover_node(int v)
    {
        m_Q.enqueue(v);
        m_discovered[v] = true;        
    }
    
                
private:
    const AdjList& m_adjlist;
    Queue<int> m_Q;
    
    Array<bool> m_discovered;
    
}; // end class BFSIterator
   


/**
 * The extended iterator of graph nodes using BFS order
 *
 * This iterator also keep tracks of predecessors and distances
 * from the seeds
 */     
class BFSIteratorEx
{
public:
        
    BFSIteratorEx(const AdjList& adjlist) 
    : m_adjlist(adjlist), m_Q(adjlist.nnodes())
    , m_discovered(adjlist.nnodes())
    , m_preds(adjlist.nnodes())
    , m_dists(adjlist.nnodes())
    {
        m_discovered.set_zeros();        
    }   
    
    
    /**
     * add a new node as seed, and tag v as discovered
     *
     * @param v the index of the node to be added as seed
     * @return whether v is added (true if v is undiscovered before adding) 
     *
     * @remarks      
     *  - this function adds v only when v remains undiscovered  
     */
    bool add_seed(int v)
    {
        if (!is_discovered(v))
        {
            discover_seed(v);
            return true;
        } 
        else
        {
            return false;
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
            m_Q.dequeue();
            
            int c = m_adjlist.neighbor_num(s);
            if (c > 0)
            {
                const int *ts = m_adjlist.neighbor_nodes(s);
                for (int i = 0; i < c; ++i)
                {
                    int t = ts[i];
                    if (!is_discovered(t))
                    {
                        discover_node(s, t);
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
    
    bool is_discovered(int v) const
    {
        return m_discovered[v];
    }
    
    int predecessor_of(int v) const
    {
        return m_preds[v];
    }
    
    int distance_of(int v) const
    {
        return m_dists[v];
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
    
    void discover_seed(int v)
    {
        m_Q.enqueue(v);
        m_discovered[v] = true;
        
        m_preds[v] = -1;
        m_dists[v] = 0;
    }
    
    void discover_node(int pred, int v)
    {
        m_Q.enqueue(v);
        m_discovered[v] = true; 
        
        m_preds[v] = pred;
        m_dists[v] = m_dists[pred] + 1;        
    }
    
                
private:
    const AdjList& m_adjlist;
    Queue<int> m_Q;
    
    Array<bool> m_discovered;
    
    Array<int> m_preds;
    Array<int> m_dists;
    
}; // end class BFSIteratorEx


enum DFSEdgeType
{
    DFS_NO_EDGE = 0,
    DFS_FORWARD_EDGE = 1,
    DFS_BACK_EDGE = 2,
    DFS_CROSS_EDGE = 3
};


struct DFSEntry
{
    int s;              // source node index
    int nnb;            // # neighbors
    const int *nbs;     // neighbors
    int i;              // offset of next neighbor

    bool remain() const
    {
        return i < nnb;
    }

    int next_neighbor()
    {
        return nbs[i++];
    }
    
    DFSEntry() { }

    DFSEntry(const AdjList& adjlist, int v)
    {
        s = v;
        nnb = adjlist.neighbor_num(v);
        nbs = adjlist.neighbor_nodes(v);
        i = 0;        
    }
};



/**
 * The iterator of graph nodes using DFS order
 */     
class DFSIterator
{    
public:
        
    DFSIterator(const AdjList& adjlist) 
    : m_adjlist(adjlist), m_seed(-1), m_stack(adjlist.nnodes())
    , m_discovered(adjlist.nnodes()), m_finished(adjlist.nnodes())
    {
        m_discovered.set_zeros();
        m_finished.set_zeros();
    }   
    
    
    /**
     * set a node as seed.
     *
     * This function acts only when v is undiscovered.  
     */
    bool set_seed(int v)
    {
        if (!is_discovered(v))
        {
            m_seed = v;
            return true;
        } 
        else
        {
            return false;
        }
    }
    
    
    /**
     * Get the next node and move the iterator forward
     *
     * @return the index of next node, or -1 when no next node
     */
    int next()
    {
        DFSEdgeType ety;
        int v = -1;
        while ((ety = next_edge(v)) && ety != DFS_FORWARD_EDGE);                
        
        return v;
    }
    
    
    /**
     * Get the next edge and move the iteratro forward
     *
     * @param t the target node of next edge (-1 if no next edge)
     *
     * @return the type of next edge
     */
    DFSEdgeType next_edge(int& t)
    {        
        while (!m_stack.empty())
        {
            DFSEntry& e = m_stack.top();
                        
            if (e.remain())
            {                                
                t = e.next_neighbor();
                                
                if (is_discovered(t))
                {
                    return is_finished(t) ? DFS_CROSS_EDGE : DFS_BACK_EDGE;
                }
                else
                {
                    discover_node(t);
                    return DFS_FORWARD_EDGE;
                }
            }
            else
            {
                finish_top();
            }
        }
        
        if (m_seed >= 0)  // starting
        {
            discover_node(m_seed);
            t = m_seed;
            m_seed = -1;
            return DFS_FORWARD_EDGE;
        }
        else  // ending
        {
            t = -1;
            return DFS_NO_EDGE;
        }
        
    }
    
    
    
public:
    bool is_discovered(int v) const
    {
        return m_discovered[v];
    }
    
    bool is_finished(int v) const
    {
        return m_finished[v];
    }
    
    int top_node() const
    {
        return m_stack.top().s;
    }
    
    
private:
    void discover_node(int v)
    {
        m_discovered[v] = true;
        m_stack.push(DFSEntry(m_adjlist, v));
    }
    
    void finish_top()
    {
        m_finished[top_node()] = true;                
        m_stack.pop();        
    }
    
                
private:
    const AdjList& m_adjlist;
    int m_seed;
    Stack<DFSEntry> m_stack;
    
    Array<bool> m_discovered;
    Array<bool> m_finished;
    
    
}; // end class DFSIterator



/**
 * The extended iterator of graph nodes using DFS order
 *
 * This class also outputs predecessors, the time-stamps
 * of discovery and finishing, as well as a sequence of nodes
 * in finishing order
 */     
class DFSIteratorEx
{    
public:
        
    DFSIteratorEx(const AdjList& adjlist) 
    : m_adjlist(adjlist), m_seed(-1), m_stack(adjlist.nnodes())
    , m_discovered(adjlist.nnodes()), m_finished(adjlist.nnodes())
    , m_time(0)
    , m_discover_time(adjlist.nnodes()), m_finish_time(adjlist.nnodes())
    , m_finish_order(adjlist.nnodes()), m_preds(adjlist.nnodes())
    {
        m_discovered.set_zeros();
        m_finished.set_zeros();
    }   
    
    
    /**
     * set a node as seed.
     *
     * This function acts only when v is undiscovered.  
     */
    bool set_seed(int v)
    {
        if (!is_discovered(v))
        {
            m_seed = v;
            return true;
        } 
        else
        {
            return false;
        }
    }
    
    
    /**
     * Get the next node and move the iterator forward
     *
     * @return the index of next node, or -1 when no next node
     */
    int next()
    {
        DFSEdgeType ety;
        int v = -1;
        while ((ety = next_edge(v)) && ety != DFS_FORWARD_EDGE);                
        
        return v;
    }
    
    
    /**
     * Get the next edge and move the iteratro forward
     *
     * @param t the target node of next edge (-1 if no next edge)
     *
     * @return the type of next edge
     */
    DFSEdgeType next_edge(int& t)
    {        
        while (!m_stack.empty())
        {
            DFSEntry& e = m_stack.top();
                        
            if (e.remain())
            {                                
                t = e.next_neighbor();
                                
                if (is_discovered(t))
                {
                    return is_finished(t) ? DFS_CROSS_EDGE : DFS_BACK_EDGE;
                }
                else
                {
                    discover_node(e.s, t);
                    return DFS_FORWARD_EDGE;
                }
            }
            else
            {
                finish_top();
            }
        }
        
        if (m_seed >= 0)  // starting
        {
            discover_node(-1, m_seed);
            t = m_seed;
            m_seed = -1;
            return DFS_FORWARD_EDGE;
        }
        else  // ending
        {
            t = -1;
            return DFS_NO_EDGE;
        }
        
    }    
    
    
public:
    bool is_discovered(int v) const
    {
        return m_discovered[v];
    }
    
    bool is_finished(int v) const
    {
        return m_finished[v];
    }
    
    int top_node() const
    {
        return m_stack.top().s;
    }
    
    int current_time() const
    {
        return m_time;
    }
    
    int discover_time_of(int i) const
    {
        return m_discover_time[i];
    }
    
    int finish_time_of(int i) const
    {
        return m_finish_time[i];
    }
    
    int predecessor_of(int i) const
    {
        return m_preds[i];
    }        
    
    const SeqList<int>& finish_order() const
    {
        return m_finish_order;
    }
    
private:
    void discover_node(int p, int v)
    {
        m_discovered[v] = true;
        m_preds[v] = p;
        
        m_discover_time[v] = next_timestamp();
        
        m_stack.push(DFSEntry(m_adjlist, v));
    }
    
    void finish_top()
    {
        int v = top_node();
        m_finished[v] = true;
        m_stack.pop();
        
        m_finish_time[v] = next_timestamp();
        m_finish_order.add(v);
    }
        
    
    int next_timestamp()
    {
       return m_time ++; 
    }
    
                
private:
    const AdjList& m_adjlist;
    int m_seed;
    Stack<DFSEntry> m_stack;
    
    Array<bool> m_discovered;
    Array<bool> m_finished;
    
    int m_time;        
    Array<int> m_discover_time;
    Array<int> m_finish_time;
    SeqList<int> m_finish_order;
    
    Array<int> m_preds;
    
}; // end class DFSIteratorEx





}


#endif

