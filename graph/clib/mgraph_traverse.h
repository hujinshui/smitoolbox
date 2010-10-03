/********************************************************************
 *  
 *  mgraph_traverse.h
 *
 *  The header files for graph traversal algorithms
 *
 *  Created by Dahua Lin, on Oct 2, 2010
 *
 ********************************************************************/

#ifndef SMI_MGRAPH_TRAVERSE_H
#define SMI_MGRAPH_TRAVERSE_H

#include "mgraph.h"
#include <string.h>

namespace smi
{
  
/**
 * trace back the path from starting vertex to a specified point
 * 
 * @param parents  the array of preceding node indices
 * @param v        the target node index
 * @param path     the buffer to output the path (in visiting order)
 *
 * @return the number of nodes on the path
 */    
int gpath_trace_back(const int *prec, int v, int *path);
    
    
/**
 * A class to support traversal algorithm and to store the 
 * traversal history
 */    
class GTraversal
{
public:
    explicit GTraversal(const GNeighborHood& G)
    : m_G(G), m_nv(0)
    {        
        int n = G.nnodes();
        m_trnodes = new int[n];
        m_visited = new bool[n];
        ::memset(m_visited, 0, n * size(bool));
        m_precs = new int[n];
    }
    
    ~GTraversal()
    {
        delete[] m_trnodes;
        delete[] m_visited;
        delete[] m_precs;
    }
    
public:
    
    const GNeighborHood& neighborhood() const
    {
        return m_G;
    }
    
    int num_nodes() const
    {
        return m_G.nnodes();
    }
    
    int num_visited() const
    {
        return m_nv;
    }
        
    bool is_visited(int v) const
    {
        return m_visited[v];
    }
    
    int prec_of(int v) const
    {
        return m_precs(v);
    }
        
    void visit(int prec, int v) 
    {
        m_visited[i] = true;
        m_trnodes[m_nv] = v;
        m_precs[m_nv] = parent;
        ++ m_nv;
    }
            
    /**
     * @param v the index of the node to trace back
     * @param path the buffer to store the visiting sequence to v
     *
     * @return the number of nodes along the path
     */
    int trace_back(int v, int *path) const
    {
        return gnode_trace_back(m_precs, v, path);
    }
    
    const int *traverse_nodes() const
    {
        return m_trnodes;
    }
    
    const bool *visited() const
    {
        return m_visited;
    }
    
    const int *preceding_nodes() const
    {
        return m_precs;
    }
                    
private:
    GTraversal(const GTraversal&);
    GTraversal& operator = (const GTraversal& );
    
private:
    const GNeighborHood& m_G;   // the base graph: (outgoing) neighborhood
    
    int m_nv;           // the number of visited nodes;
    int *m_trnodes;     // the pre-allocated buffer to store the sequence of visited nodes
    int *m_visited;     // the indicators of whether the nodes are visited.
    int *m_precs;       // the preceding nodes of visited nodes    
    
}; // end class GTraversal
    
    

/**
 * The base class for graph traversal
 */
class GTraverserBase
{
public:
    GTraverserBase(const GNeighborHood& G) : m_traversal(G)
    {
    }
    
    virtual ~GTraverserBase()
    {
    }
    
    const GTraversal& traversal() const
    {
        return m_traversal;
    }
           
    /**
     * initialize the internal data structure 
     * 
     * @param ns the number of seed nodes
     * @param seeds the array of seed node indices
     */
    virtual void initialize(int ns, int *seeds) = 0;    
            
    /**
     * Performs search along the graph to a stop node
     *
     * @param vstop stops the search when vstop is visited 
     *        (set vstop to -1 for traversing the entire graph)
     * @param amap the indicators of accessibility 
     *        (set amap to NULL when all nodes are accessible)
     *
     * @return whether vstop has been visited
     */
    virtual bool search(int vstop, const bool *amap) = 0;
    
     
    /**
     * Performs search along the graph until all nodes in vstops
     * have been visited (or the entire connected component has been
     * visited)
     *
     * @param n the number of nodes in vstops
     * @param vstops the set of nodes such that the search stops when
     *        all of them are visited
     * 
     * @return whether all nodes in vstops are accessible
     */
    bool search_all(int n, int *vstops)
    {
        bool succeed = true;
        
        for (int i = 0; i < n; ++i)
        {            
            int v = vstops[i];
            if (!m_traversal.is_visited(v) && !this->search(v))
                succeed = false;
        }
        
        return succeed;
    }
    
    
    /**
     * Performs the traversal over the entire connected component
     */
    void traverse(const bool *amap = 0)
    {
        search(-1, amap);
    }
    
protected:
    GTraversal m_traversal;
    
}; // end class GraphSearchBase


/**
 * The class for breadth first search
 */
class BreadthFirstSearch : GraphSearchBase
{
public:
    BreadthFirstSearch(const GNeighborHood& G) 
    : GraphSearchBase(G), _impl(0) { }
    
    ~BreadthFirstSearch()
    {
        if (_impl != 0) delete _impl;
    }
                
    virtual void initialize(int ns, int *starts);
              
    virtual bool search(int vstop, const bool *amap);
            
private:
    
    class Impl;
    Impl* _impl;
};



/**
 * The class for depth first search
 */
class DepthFirstSearch : GraphSearchBase
{
public:
    DepthFirstSearch(const GNeighborHood& G) 
    : GraphSearchBase(G), _impl(0) { }
    
    ~DepthFirstSearch()
    {
        if (_impl != 0) delete _impl;
    }
                
    virtual void initialize(int ns, int *starts);    
    
    virtual bool search(int vstop, const bool *amap);
            
private:
    
    class Impl;
    Impl* _impl;
};

}


#endif

