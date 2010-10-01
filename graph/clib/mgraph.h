/********************************************************************
 *
 *  mgraph.h
 *
 *  The class to provide a lightweight C++ representation
 *  of graph.
 *
 *  Created by Dahua Lin, on Oct 1, 2010
 *
 *******************************************************************/

#ifndef SMI_MGRAPH_H
#define SMI_MGRAPH_H

#include <mex.h>

namespace smi
{
    
struct Edge
{
    int s;  // source index
    int t;  // target index
    
    Edge() { }
    
    Edge(int _s, int _t): s(_s), t(_t) { }    
};


template<typename T>
struct WEdge
{
    int s;  // source index
    int t;  // target index
    T w;    // the weight
    
    WEdge() { }
    
    WEdge(int _s, int _t, T _w) : s(_s), t(_t), w(_w) { }
};


template<typename T>
class GNeighborHood
{
public:
    typedef T weight_type;
    
public:
    GNeighborHood(int n, int m, const int *I, const int *J, const weight_type *W)
    {
        // construct the neighborhood structure
        
        m_nnodes = n;
        m_nedges = m;
        
        // prepare storage
        
        m_offsets = new int[n];
        m_nnbs = new int[n];
                
        m_J = new int[m];
        m_W = new weight_type[m];
        
        // count neighbors
        
        for (int i = 0; i < n; ++i)
        {
            m_nnbs[k] = 0;
        }
        
        for (int k = 0; k < m; ++k)
        {
            ++ m_nnbs[I[k]];
        }
        
        // calculate offsets
        
        int *c = new int[n];
        
        int p = 0;
        for (int i = 0; i < n; ++i, p += m_nnbs[i])
        {
            m_offsets[i] = c[i] = p;            
        }
        
        // fill in neighbors
        
        for (int k = 0; k < m; ++k)
        {
            int j = c[I[k]]++;            
            m_J[j] = J[k]; 
            m_W[j] = W[k];
        }        
        
        delete[] c;                
    }
    
    ~GNeighborHood()
    {
        delete[] m_offsets;
        delete[] m_nnbs;
        delete[] m_J;
        delete[] m_W;
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
    
    const int* neighbor_indices(int i) const
    {
        return m_J + m_offsets[i];
    }
    
    const int* neighbor_weights(int i) const
    {
        return m_W + m_offsets[i];
    }
    
    int neighbor_num(int i) const
    {
        return m_nnbs[i];
    }
    
private:
    int m_nnodes;    // the number of nodes
    int m_nedges;    // the number of edges
    int *m_offsets;  // offset for each source: length n
    int *m_nnbs;     // # neighbors for each source: length n
    int *m_J;        // concatenated array of neighbor indices: length m
    int *m_W;        // concatenated array of neighbor weights: length m
    
}; // end class GNeighborHood


enum GNeighborHoodType
{
    GNB_OUT,
    GNB_INT
};


// a light weight interface for matlab input of graph    
template<typename T>
class MGraph
{
public:
    typedef T weight_type;
    
public:
    
    MGraph(int n, int m, const int *I, const int *J, const weight_type *W)
    {
        m_nnodes = n;
        m_nedges = m;
        m_own = false;
        
        m_I = I;
        m_J = J;
        m_W = W;
    }
    
    MGraph(const mxArray *mxG)
    {
        const mxArray *mxN = mxGetField(mxG, 0, "n");
        const mxArray *mxI = mxGetField(mxG, 0, "I");
        const mxArray *mxJ = mxGetField(mxG, 0, "J");
        const mxArray *mxW = mxGetField(mxG, 0, "W");                
        
        m_nnodes = n;
        int m = mxGetNumberOfElements(mxI);
        m_nedges = m;
        m_own = false;
        
        m_I = (const int*)mxGetData(mxI);
        m_J = (const int*)mxGetData(mxJ);
        m_W = (const weight_type*)mxGetData(mxW);
    }        
    
    
    ~MGraph()
    {        
        if (m_own)
        {
            delete[] m_I;     
            delete[] m_J;
            delete[] m_W;
        }
    }
    
        
    int nnodes() const
    {
        return m_nnodes;
    }
    
    int nedges() const
    {
        return m_nedges;
    }
    
    bool own_edge_memory() const
    {
        return m_own;
    }
    
    Edge edge(int i) const
    {
        return Edge(m_I[i], m_J[i]);
    }
    
    WEdge<weight_type> wedge(int i) const
    {
        return WEdge<weight_type>(m_I[i], m_J[i], m_W[i]);
    }
    
    GNeighborHood* create_neighborhood(GNeighborHoodType nbType) const
    {
        if (nbType == GNB_OUT)
        {
            return new GNeighborHood(m_nnodes, m_nedges, m_I, m_J, m_W);
        }
        else if (nbType == GNB_IN)
        {
            return new GNeighborHood(m_nnodes, m_nedges, m_J, m_I, m_W);
        }
        else
        {
            return 0;
        }
    }
        
    GNeighborHood* create_in_neighborhood() const
    {
        return new GNeighborHood(m_nnodes, m_nedges, m_J, m_I, m_W);        
    }
                            
private:
    MGraph(const MGraph& );
    MGraph& operator = (const MGraph& );    
        
    
private:
    int m_nnodes;   // the number of nodes
    int m_nedges;   // the number of edges
        
    bool m_own;     // whether it owns the storage of I, J, and W.
    const int *m_I;         // the array of source indices (zero-based)
    const int *m_J;         // the array of dest indices (zero-based)
    const weight_type *m_W; // the array of weights    

}; // end class MGraph
    
            
}





