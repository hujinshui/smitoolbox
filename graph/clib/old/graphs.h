/********************************************************************
 *
 *  graphs.h
 *
 *  The classes to represent Graphs
 *
 *  Created by Dahua Lin, on Oct 3, 2010
 *
 ********************************************************************/

#ifndef SMI_GRAPHS_H
#define SMI_GRAPHS_H

#include <string.h>


namespace smi
{
 
struct Edge
{
    int s;  // the source node
    int t;  // the target node
    
    Edge() { }
    
    Edge(int s_, int t_) : s(s_), t(t_) { }
};


template<typename TWeight>
struct WEdge : public Edge
{
    TWeight w;  // the edge weight
    
    WEdge() { }
    
    WEdge(int s_, int t_, TWeight w_) : Edge(s_, t_), w(w_) { }
};
    
    
/**
 * The class to represent a 
 * which refers to some external memory for I and J
 */
class RefGraph
{
public:
    RefGraph(int n, int m, const int *I, const int *J)
    : m_nnodes(n), m_nedges(m), m_I(I), m_J(J) { }    
    
    int nnodes() const
    {
        return m_nnodes;
    }
    
    int nedges() const
    {
        return m_nedges;
    }
    
    
    Edge edge(int k) const
    {
        return Edge(m_I[k], m_J[k]);
    }    
    
    
    int source_of(int k) const
    {
        return m_I[k];
    }
    
    int target_of(int k) const
    {
        return m_J[k];
    }
    
    
    const int *edge_sources() const
    {
        return m_I;
    }
    
    const int *edge_targets() const
    {
        return m_J;
    }    
    
       
protected:
    int m_nnodes;
    int m_nedges;
    
    const int *m_I;
    const int *m_J;    
    
}; // end class Graph


template<typename TWeight>
class RefWGraph : public RefGraph
{
public:
    RefWGraph(int n, int m, const int *I, const int *J, const TWeight *W)
    : RefGraph(n, m, I, J), m_W(W) { }
    
    
    TWeight weight_of(int k) const
    {
        return m_W[k];
    }
    
    WEdge<TWeight> wedge(int k) const
    {
        return WEdge<TWeight>(this->m_I[k], this->m_J[k], m_W[k]);
    }    
    
    const TWeight *edge_weights() const
    {
        return m_W;
    }
    
    
private:
    const TWeight *m_W;
    
}; // end class RefWGraph


inline RefGraph transpose(const RefGraph& G)
{
    return RefGraph(G.nnodes(), G.nedges(), G.edge_targets(), G.edge_sources());
}

template<typename TWeight>
inline RefWGraph<TWeight> transpose(const RefWGraph<TWeight>& G)
{
    return RefWGraph<TWeight>(G.nnodes(), G.nedges(), 
            G.edge_targets(), G.edge_sources(), G.edge_weights());
}




    
/**
 * The class to represent an adjacency List
 */    
class AdjList
{    
public:
    AdjList(int n, int m, const int *I, const int* J)    
    {
        _init_IJ(n, m, I, J);
    }
    
    AdjList(const RefGraph& G)
    {
        _init_IJ(G.nnodes(), G.nedges(), G.edge_sources(), G.edge_targets());
    }
    
    virtual ~AdjList()
    {
        delete[] m_nnbs;
        delete[] m_offsets;
        delete[] m_nbnodes;
    }
    
    
public:
    
    int nnodes() const
    {
        return m_nnodes;
    }
    
    int nedges() const
    {
        return m_nedges;
    }
    
    int neighbor_num(int v) const
    {
        return m_nnbs[v];
    }
    
    const int *neighbor_nodes(int v) const
    {
        return m_nbnodes + m_offsets[v];
    }
    
        
protected:
    
    AdjList() { }    
    
    
    void _init_counts_and_offsets(int n, int m, const int *I)
    {                
        // construct the neighborhood structure
        
        m_nnodes = n;
        m_nedges = m;
        
        // prepare storage
        
        m_nnbs = new int[n];
        m_offsets = new int[n];        
        
        // count neighbors
        
        for (int i = 0; i < n; ++i) 
        {
            m_nnbs[i] = 0;
        }
        
        for (int k = 0; k < m; ++k) 
        {
            ++ m_nnbs[I[k]];
        }
        
        // calculate offsets
                
        int p = 0;
        for (int i = 0; i < n; ++i) 
        {
            m_offsets[i] = p;
            p += m_nnbs[i];
        }
    }
    
    
    void _init_IJ(int n, int m, const int *I, const int *J)
    {
        _init_counts_and_offsets(n, m, I);
        
        // fill in neighbors
        
        int *c = new int[n];
        ::memcpy(c, m_offsets, n * sizeof(int));        
        
        m_nbnodes = new int[m];
        
        for (int k = 0; k < m; ++k) 
        {
            int i = I[k];                        
            m_nbnodes[c[i]++] = J[k];
        }
        
        delete[] c;
    }
        
    
protected:    
    int m_nnodes;   // the number of nodes (n)
    int m_nedges;   // the number of edges (m)
    
    int *m_nnbs;    // the number of neighbors of each node (length n)
    int *m_offsets; // the offsets at m_nbs of each node (length n)        
    int *m_nbnodes; // the concatenated array of neighboring nodes (length m)    
    
}; // end class AdjList
    


template<typename TWeight>
class WAdjList : public AdjList
{
public:
    WAdjList(int n, int m, const int *I, const int *J, const TWeight *W)
    {
        _init_IJW(n, m, I, J, W);
    }
    
    WAdjList(const RefWGraph<TWeight>& G)
    {
        _init_IJW(G.nnodes(), G.nedges(), G.edge_sources(), G.edge_targets(), G.edge_weights());
    }
    
    
    virtual ~WAdjList()
    {
        delete[] m_nbweights;
    }
    
    const TWeight *neighbor_weights(int v) const
    {
        return m_nbweights + this->m_offsets[v];
    }
    
private:
    
    void _init_IJW(int n, int m, const int *I, const int *J, const TWeight *W)
    {
        this->_init_counts_and_offsets(n, m, I);
        
        // fill in neighbors
        
        int *c = new int[n];
        ::memcpy(c, this->m_offsets, n * sizeof(int));        
        
        this->m_nbnodes = new int[m];
        this->m_nbweights = new TWeight[m];
        
        for (int k = 0; k < m; ++k) 
        {
            int i = I[k];                        
            this->m_nbnodes[c[i]] = J[k];
            this->m_nbweights[c[i]++] = W[k];
        }
        
        delete[] c;        
    }
    
    
private:
    TWeight *m_nbweights;   // the concatenated array of neighboring weights (length m)
    
}; // end class WAdjList
        

}


#endif


