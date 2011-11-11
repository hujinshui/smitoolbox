/********************************************************************
 *
 *  smi_graph.h
 *
 *  A simple graph class for SMI toolbox, which has convenient
 *  interface to interact with MATLAB.
 *
 *  Note: the node/edge index uses one-based system.
 *
 *  Created by Dahua Lin, on Oct 27, 2011
 *
 ********************************************************************/

#include <stdint.h>

#include <exception>
#include <stdexcept>


namespace smi
{
    
typedef int32_t gindex_t;    

struct DEdge
{
    gindex_t s;
    gindex_t t;
};

inline DEdge make_dedge(gindex_t s, gindex_t t)
{
    DEdge ed;
    ed.s = s;
    ed.t = t;
    return ed;
}
    

struct DGraphSpec
{
    gindex_t n;
    gindex_t m;
    const DEdge *edges;
    
    const gindex_t *out_nbs;
    const gindex_t *out_edges;
    const gindex_t *out_degs;
    const gindex_t *offsets;
};


/**
 * The class to represent a directed graph
 */
class DGraph
{
public:
    DGraph(gindex_t n, gindex_t max_m)
    : m_nnodes(n)
    , m_nedges(0)
    , m_max_nedges(max_m)
    , m_own(true)
    , m_has_nbs(false)
    {
        m_edgelist = new DEdge[max_m];
        
        m_out_neighbors = 0;
        m_out_edges = 0;
        m_out_degs = 0;
        m_offsets = 0;
    }
    
    DGraph(const DGraphSpec& s)  // for MATLAB only
    : m_nnodes(s.n)
    , m_nedges(s.m)
    , m_max_nedges(s.m)
    , m_own(false)
    , m_has_nbs(s.out_nbs != 0)
    , m_edgelist(const_cast<DEdge*>(s.edges))
    , m_out_neighbors(const_cast<gindex_t*>(s.out_nbs))
    , m_out_edges(const_cast<gindex_t*>(s.out_edges))
    , m_out_degs(const_cast<gindex_t*>(s.out_degs))
    , m_offsets(const_cast<gindex_t*>(s.offsets))
    {
    }
        
    
    ~DGraph()
    {
        if (m_own)
        {
            if (m_edgelist) delete [] m_edgelist;
            if (m_out_neighbors) delete [] m_out_neighbors;
            if (m_out_edges) delete [] m_out_edges;
            if (m_out_degs) delete [] m_out_degs;
            if (m_offsets) delete [] m_offsets;
        }
    }
    
private:
    // disable copying
    
    DGraph(const DGraph& );    
    DGraph& operator = (const DGraph& );
    
public:
    // basic information
    
    gindex_t nnodes() const
    {
        return m_nnodes;
    }
    
    gindex_t nedges() const
    {
        return m_nedges;
    }
    
    gindex_t max_nedges() const
    {
        return m_max_nedges;
    }
    
    const DEdge& edge(gindex_t i) const  // one-based index
    {
        return m_edgelist[i-1];
    }
    
    bool has_neighborhood_system() const
    {
        return m_has_nbs;
    }
    
    DGraphSpec get_spec() const
    {
        DGraphSpec s;
        
        s.n = m_nnodes;
        s.m = m_nedges;
        s.edges = m_edgelist;
        
        s.out_nbs = m_out_neighbors;
        s.out_edges = m_out_edges;
        s.out_degs = m_out_degs;
        s.offsets = m_offsets;
        
        return s;
    }
        
public:
    
    // NOTE: the following function requires neighborhoods.
    // Please make sure the neighborhood is ready before using them,
    // otherwise it may lead to undefined behavior
    
    gindex_t out_degree(gindex_t v) const
    {
        return m_out_degs[v-1];
    }
    
    const gindex_t* out_neighbors(gindex_t v) const
    {
        return m_out_neighbors + m_offsets[v-1];
    }
    
    const gindex_t* out_edges(gindex_t v) const
    {
        return m_out_edges + m_offsets[v-1];
    }
        
    
public:
    // building
    
    void add_edge(gindex_t s, gindex_t t)
    {
        add_edge(make_dedge(s, t));
    }
    
    void add_edge(const DEdge& e) 
    {
        if (m_has_nbs)
        {
            throw std::runtime_error(
                    "Cannot modify the graph when neighborhood has been established.");
        }
        
        if (m_nedges == m_max_nedges)
        {
            throw std::runtime_error(
                    "The edge list is full.");
        }
        
        m_edgelist[m_nedges++] = e;        
    }
    
    void build_neighborhood_system()
    {
        if (m_has_nbs)
        {
            throw std::runtime_error(
                    "The neighborhood has been established.");
        }
        
        gindex_t n = m_nnodes;
        gindex_t m = m_nedges;
        
        // pass 1: set up degrees and offsets
        
        m_out_degs = new gindex_t[n];
        m_offsets  = new gindex_t[n];
        
        for (gindex_t i = 0; i < n; ++i)
        {
            m_out_degs[i] = 0;
        }
        
        for (gindex_t i = 0; i < m; ++i)
        {
            ++ m_out_degs[m_edgelist[i].s - 1]; // edge has one-based index
        }
        
        gindex_t co = 0;
        for (gindex_t i = 0; i < n; ++i)
        {
            m_offsets[i] = co;
            co += m_out_degs[i];
        }
        
        // pass 2: fill in neighbors and incident edges
        
        m_out_neighbors = new gindex_t[m];
        m_out_edges = new gindex_t[m];
        
        gindex_t *cs = new gindex_t[n];
        for (gindex_t i = 0; i < n; ++i)
        {
            cs[i] = 0;
        }        
        
        for (gindex_t i = 0; i < m; ++i)
        {             
            const DEdge& e = m_edgelist[i];           
            gindex_t cc = cs[e.s-1] ++;
            int j = m_offsets[e.s-1] + cc;
            
            m_out_neighbors[j] = e.t;
            m_out_edges[j] = i + 1;  // one-based index
        }
        
        m_has_nbs = true;
        
    }
            
    
private:
    gindex_t m_nnodes;
    gindex_t m_nedges;
    gindex_t m_max_nedges;
    
    bool m_own;  // whether it owns edge memory
    bool m_has_nbs; // whether the neighborhood system has been built
    
    DEdge *m_edgelist;            
    gindex_t *m_out_neighbors;  // length = m
    gindex_t *m_out_edges;      // length = m
    gindex_t *m_out_degs;       // length = n
    gindex_t *m_offsets;        // length = n
    
}; // end class Graph
    
        
}


