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

#include <stack>
#include <queue>
#include <valarray>

namespace smi
{
  
/**
 * The iterator of graph nodes using BFS order
 */     
class BFSIterator
{    
public:
    struct entry
    {
        int pred;   // the predecessor of node
        int node;   // the currently visiting node
        
        entry() { }        
        entry(int p, int v) : pred(p), node(v) { }
    };
    
    
    struct nil_observer_t
    {
        void on_discover(int p, int v) { }
        void on_visit (int p, int v) { }
        
    } nil_observer;
    
public:            
    BFSIterator(const AdjList& G) 
    : m_adjlist(G), m_Q(), m_discovered(G.nnodes())
    {
        m_discovered.set_zeros();
    }   
    
    
    int num_nodes() const
    {
        return m_adjlist.nnodes();
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
            discover_node(-1, v);
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
     * @param obs the observer that reacts on the next step
     *
     * @returns the index of next visiting node, -1 if no next node.
     */
    template<typename TObserver>
    int next(TObserver& obs)
    {
        if (!m_Q.empty())
        {
            // visit the next node
            const entry& e = m_Q.front();            
            obs.on_visit(e.pred, e.node);                        
            int s = e.node;
            m_Q.pop();
            
            // add its unvisited neighbors to queue            
            int c = m_adjlist.neighbor_num(s);
            if (c > 0)
            {
                const int *ts = m_adjlist.neighbor_nodes(s);
                for (int i = 0; i < c; ++i)
                {
                    int t = ts[i];
                    if (!is_discovered(t))
                    {
                        obs.on_discover(s, t);
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
    
    
    int next()
    {
        return next(nil_observer);
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
    
    void discover_node(int pred, int v)
    {
        m_Q.push(entry(pred, v));
        m_discovered[v] = true;                    
    }        
    
                
private:
    const AdjList& m_adjlist;
    std::queue<entry> m_Q;    
    BoolArray m_discovered;
    
}; // end class BFSIterator
   



enum DFSEdgeType
{
    DFS_NO_EDGE = 0,
    DFS_FORWARD_EDGE = 1,
    DFS_BACK_EDGE = 2,
    DFS_CROSS_EDGE = 3
};


/**
 * The iterator of graph nodes using DFS order
 */     
class DFSIterator
{    
public:
    
    struct edge_info
    {
        int s;
        int t;
        DFSEdgeType etype;         
        
        edge_info() { }
        edge_info(int s_, int t_, DFSEdgeType ety) : s(s_), t(t_), etype(ety) { }        
    };
    
    
    struct entry
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

        entry() { }

        entry(const AdjList& adjlist, int v)
        {
            s = v;
            nnb = adjlist.neighbor_num(v);
            nbs = adjlist.neighbor_nodes(v);
            i = 0;        
        }
    };
    
    
    struct nil_observer_t
    {
        void on_discover(int s, int t) { }
        void on_finish(int v) { }
        
    } nil_observer;
    
    
public:
        
    DFSIterator(const AdjList& G) 
    : m_adjlist(G), m_seed(-1), m_stack()
    , m_discovered(G.nnodes()), m_finished(G.nnodes())
    {
        m_discovered.set_zeros();
        m_finished.set_zeros();
    }   
    
    
    int num_nodes() const
    {
        return m_adjlist.nnodes();
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
     * @param the observer that reacts to DFS steps
     *
     * @return the index of next node, or -1 when no next node
     */
    template<typename TObserver>
    int next(TObserver& obs)
    {        
        edge_info e = next_edge(obs);
        while (e.etype != DFS_FORWARD_EDGE && e.etype != DFS_NO_EDGE) 
            e = next_edge(obs);
                        
        return e.t;
    }
    
    
    int next()
    {
        return next(nil_observer);
    }
    
    
    
    /**
     * Get the next edge and move the iteratro forward
     * 
     * @return the information struct of next edge
     */
    template<typename TObserver>
    edge_info next_edge(TObserver& obs)
    {        
        while (!m_stack.empty())
        {
            entry& e = m_stack.top();
                        
            if (e.remain())
            {                                
                int t = e.next_neighbor();
                                
                if (is_discovered(t))
                {
                    DFSEdgeType ety = is_finished(t) ? DFS_CROSS_EDGE : DFS_BACK_EDGE;                            
                    return edge_info(e.s, t, ety);                        
                }
                else
                {
                    obs.on_discover(e.s, t);
                    discover_node(t);
                    
                    return edge_info(e.s, t, DFS_FORWARD_EDGE);
                }
            }
            else
            {
                int vtop = top_node();
                
                obs.on_finish(vtop);
                m_finished[vtop] = true;                
                m_stack.pop();  
            }
        }
        
        if (m_seed >= 0)  // starting
        {
            obs.on_discover(-1, m_seed);
            discover_node(m_seed);
            
            int t = m_seed;
            m_seed = -1;
            return edge_info(-1, t, DFS_FORWARD_EDGE);
        }
        else  // ending
        {
            return edge_info(-1, -1, DFS_NO_EDGE);
        }
        
    }
    
    
    edge_info next_edge()
    {
        return next_edge(nil_observer);
    }    
            
    
    /**
     * Move the iterator forward until a loop is detected (back edge)
     *
     * @return the target node of the back edge, or -1 if no loop is found
     */
    template<typename TObserver>
    int detect_loop(TObserver& obs)
    {
        edge_info e = next_edge(obs);
        while (e.etype != DFS_BACK_EDGE && e.etype != DFS_NO_EDGE) 
            e = next_edge(obs);
                        
        return e.t;
    }
    
    int detect_loop()
    {
        detect_loop(nil_observer);
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
        m_stack.push(entry(m_adjlist, v));
    }        
                
private:
    const AdjList& m_adjlist;
    int m_seed;
    std::stack<entry> m_stack;
    
    BoolArray m_discovered;
    BoolArray m_finished;
    
    
}; // end class DFSIterator




/*****************************************
 *
 *  Some algorithms based on BFS or DFS
 *
 *****************************************/

template<typename TInserter>
struct acyclic_testing_dfs_observer
{
    acyclic_testing_dfs_observer(TInserter it) : inserter(it) { }
    
    void on_discover(int s, int t) { }
    
    void on_finish(int v)
    {
        *(inserter++) = v;
    }
    
    TInserter inserter;
};



/**
 * tests whether a directed graph is acylic
 *
 * @param dfs_it A DFS iterator (at initial status) on the graph 
 * @param inserter for output the sequence of finished node
 * 
 * @return true if the graph is acyclic, false otherwise
 *
 * @remarks if the graph is acyclic, the reversed order of the 
 *  finishing order is a topological order of the graph
 *
 */
template<typename TInserter>
bool test_acyclic(DFSIterator& dfs_it, TInserter inserter)
{           
    int n = dfs_it.num_nodes();
    acyclic_testing_dfs_observer<TInserter> obs(inserter);
    
    for (int i = 0; i < n; ++i)
    {
        if (!dfs_it.is_discovered(i))
        {
            dfs_it.set_seed(i);
            
            if (dfs_it.detect_loop(obs) >= 0)
                return false;
        }
    }
    return true;
}


/**
 * Get connected components of an undirected (or symmetric) graph
 *
 * @param adjlist the adjacency list of the graph (input)
 * @param it_csizs the output iterator of sizes (# nodes) of components
 * @param it_nodes the output iterator of nodes in each component
 *
 * @return the number of components
 *
 * @remarks it uses BFS to traverse each component
 */
template<typename TIterSize, typename TIterNodes>
int get_connected_components(const AdjList& adjlist, 
        TIterSize it_csizs, TIterNodes it_nodes)
{
    int n = adjlist.nnodes();
    
    BFSIterator bfs_it(adjlist);
        
    int c = 0;        
    for (int i = 0; i < n; ++i)
    {
        if (!bfs_it.is_discovered(i))
        {
            bfs_it.add_seed(i);
            
            int v;
            int cs = 0;
            while ((v = bfs_it.next()) >= 0)
            {
                *(it_nodes++) = v; 
                ++ cs;
            }
        
            *(it_csizs++) = cs;
            ++c;
        }
    }
    
    return c;
}


}


#endif

