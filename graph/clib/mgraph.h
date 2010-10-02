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


/*************************************
 *
 *  Graph classes
 *
 ************************************/

 
class MGraph
{    
public:    
    static mxClassID get_weight_class(const mxArray *mxG)
    {
        const mxArray *mxW = mxGetField(mxG, 0, "W");
        return mxGetClassID(mxW);
    }
    
public:
    
    MGraph(int n, int m, const int *I, const int *J)
    {
        m_nnodes = n;
        m_nedges = m;
        
        m_own = false;        
        m_I = I;
        m_J = J;
    }
    
    MGraph(const mxArray *mxG)
    {
        const mxArray *mxN = mxGetField(mxG, 0, "n");
        const mxArray *mxI = mxGetField(mxG, 0, "I");
        const mxArray *mxJ = mxGetField(mxG, 0, "J");            
        
        m_nnodes = (int)mxGetScalar(mxN);
        int m = mxGetNumberOfElements(mxI);
        m_nedges = m;
        
        m_own = false;        
        m_I = (const int*)mxGetData(mxI);
        m_J = (const int*)mxGetData(mxJ);
    }        
    
    
    virtual ~MGraph()
    {        
        if (m_own)
        {
            delete[] m_I;     
            delete[] m_J;
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
    
public:
    const int *edge_sources() const
    {
        return m_I;
    }
    
    const int *edge_targets() const
    {
        return m_J;
    }
                                                    
private:
    MGraph(const MGraph& );
    MGraph& operator = (const MGraph& );            
    
protected:
    int m_nnodes;   // the number of nodes
    int m_nedges;   // the number of edges
        
    bool m_own;     // whether it owns the storage of I, J, and W.
    const int *m_I;         // the array of source indices (zero-based)
    const int *m_J;         // the array of dest indices (zero-based) 

}; // end class MGraph


template<typename T>
class MWGraph : public MGraph
{
public:
    typedef T weight_type;
    
    MWGraph(int n, int m, const int *I, const int *J, const weight_type *W)
    : MGraph(n, m, I, J)
    {
        m_W = W;
    }
    
    MWGraph(const mxArray *mxG) : MGraph(mxG)
    {
        const mxArray *mxW = mxGetField(mxG, 0, "W");
        m_W = (const weight_type*)mxGetData(mxW);
    }
    
    virtual ~MWGraph()
    {
        if (this->m_own)
        {
            delete[] m_W;
        }
    }
    
public:
    WEdge<weight_type> wedge(int i) const
    {
        return WEdge<weight_type>(this->m_I[i], this->m_J[i], this->m_W[i]);
    }
    
    const weight_type *edge_weights() const
    {
        return m_W;
    }
                
private:
    const weight_type *m_W;  // the array of weights
    
}; // end class MWGraph
  


/*************************************
 *
 *  Neighborhood system
 *
 ************************************/

struct gnb_out { };
struct gnb_in { };


class GNeighborHood
{
public:
    GNeighborHood(int n, int m, const int *I, const int *J)
    {
        _init(n, m, I, J);
    }
    
    GNeighborHood(const MGraph& G0, gnb_out)
    {
        _init(G0.nnodes(), G0.nedges(), G0.edge_sources(), G0.edge_targets());
    }
    
    GNeighborHood(const MGraph& G0, gnb_in)
    {
        _init(G0.nnodes(), G0.nedges(), G0.edge_targets(), G0.edge_sources());
    }
    
    virtual ~GNeighborHood()
    {
        delete[] m_offsets;
        delete[] m_nnbs;
        delete[] m_J;
        delete[] m_E;
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
    
    const int* neighbor_nodes(int i) const
    {
        return m_J + m_offsets[i];
    }
    
    const int* neighbor_edges(int i) const
    {
        return m_E + m_offsets[i];
    }
        
    int neighbor_num(int i) const
    {
        return m_nnbs[i];
    }
    
    int max_neighbor_num() const
    {
        return m_max_nnbs;
    }
    
public:
    virtual void mdump(int ibase=1)
    {
        mexPrintf("# nodes = %d\n", this->nnodes());
        mexPrintf("# edges = %d\n", this->nedges());
        mexPrintf("------------------------\n");
        for (int i = 0; i < m_nnodes; ++i)
        {
            mexPrintf("[%d]: ", ibase + i);
            
            int nnb = this->neighbor_num(i);
            const int *nbs = this->neighbor_nodes(i);
            
            for (int j = 0; j < nnb; ++j)
            {
                mexPrintf("%d ", ibase + nbs[j]);
            }
            mexPrintf("\n");
        }        
    }
    
private:
    GNeighborHood(const GNeighborHood& );
    GNeighborHood& operator = (const GNeighborHood& );   
    
private:
    void _init(int n, int m, const int *I, const int *J)
    {
        // construct the neighborhood structure
        
        m_nnodes = n;
        m_nedges = m;
        m_max_nnbs = 0;
        
        // prepare storage
        
        m_offsets = new int[n];
        m_nnbs = new int[n];
                
        m_J = new int[m];
        m_E = new int[m];
        
        // count neighbors
        
        for (int i = 0; i < n; ++i)
        {
            m_nnbs[i] = 0;
        }
        
        for (int k = 0; k < m; ++k)
        {
            ++ m_nnbs[I[k]];
        }
        
        for (int i = 0; i < n; ++i)
        {
            if (m_nnbs[i] > m_max_nnbs) 
                m_max_nnbs = m_nnbs[i];
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
            m_E[j] = k;
        }        
        
        delete[] c;      
    }
    
protected:
    int m_nnodes;    // the number of nodes
    int m_nedges;    // the number of edges
    int m_max_nnbs;  // maximum $nbs of a node
    
    int *m_offsets;  // offset for each source: length n
    int *m_nnbs;     // # neighbors for each source: length n
    int *m_J;        // concatenated array of neighbor indices: length m
    int *m_E;        // concatenated array of corresponding edge indices
    
}; // end class GNeighborHood


template<typename T>
class WGNeighborHood : public GNeighborHood
{
public:
    typedef T weight_type;
    
public:
    WGNeighborHood(int n, int m, const int *I, const int *J, const weight_type *W)
    : GNeighborHood(n, m, I, J)
    {
        _init_W(W);
    }
    
    WGNeighborHood(const MWGraph<T>& G0, gnb_out)
    : GNeighborHood(G0, gnb_out())
    {
        _init_W(G0.edge_weights());
    }
    
    WGNeighborHood(const MWGraph<T>& G0, gnb_in)
    : GNeighborHood(G0, gnb_in())
    {
        _init_W(G0.edge_weights());
    }
    
    virtual ~WGNeighborHood()
    {
        delete[] m_W;
    }
    
    const weight_type* neighbor_weights(int i) const
    {
        return m_W + this->m_offsets[i];
    }
    
    virtual void mdump(const char *wfmt, int ibase=1)
    {
        mexPrintf("# nodes = %d\n", this->nnodes());
        mexPrintf("# edges = %d\n", this->nedges());
        mexPrintf("------------------------\n");
        for (int i = 0; i < m_nnodes; ++i)
        {
            mexPrintf("[%d]: ", ibase + i);
            
            int nnb = this->neighbor_num(i);
            const int *nbs = this->neighbor_nodes(i);
            const weight_type *ws = this->neighbor_weights(i);
            
            for (int j = 0; j < nnb; ++j)
            {
                mexPrintf("(%d, ", ibase + nbs[j]);
                mexPrintf(wfmt, ws[j]);
                mexPrintf(")");
            }
            mexPrintf("\n");
        }        
    }
    
private:
    void _init_W(const weight_type *W)
    {
        int m = this->m_nedges;
        m_W = new weight_type[m];        
        const int *E = this->m_E;
        
        for (int k = 0; k < m; ++k)
        {
            m_W[k] = W[E[k]];
        }
    }
 
private:
    weight_type *m_W;   // concatenated array of neighbor weights: length m
    
}; // end class WGNeighborHood


}

#endif





