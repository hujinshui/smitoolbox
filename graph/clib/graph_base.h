/********************************************************************
 *
 *  graph_base.h
 *
 *  The BGL-compliant graph classes which are designed for efficient
 *  adaptation from MATLAB representation.
 *  
 *  Created by Dahua Lin, on Oct 8, 2010
 *
 ********************************************************************/

#ifndef SMI_GRAPH_BASE_H
#define SMI_GRAPH_BASE_H


#include <boost/iterator/iterator_facade.hpp>
#include <boost/graph/graph_traits.hpp>

namespace smi
{

// Basic types    
    
typedef int graph_size_t;        

struct vertex_t
{
    graph_size_t i;
    
    vertex_t { }
    vertex_t (graph_size_t i_) : i(i_) { }
    
    bool operator == (const vertex_t& rhs) const
    {
        return i == rhs.i;
    }
    
    bool operator != (const vertex_t& rhs) const
    {
        return i != rhs.i;
    }    
};


struct edge_t
{
    vertex_t s;
    vertex_t t;     
    
    edge_t() { }
    edge_t(vertex_t s_, vertex_t t_) : s(s_), t(t_) { }
    
    bool operator == (const edge_t& rhs) const
    {
        return s == rhs.s && t == rhs.t;
    }
    
    bool operator != (const edge_t& rhs) const
    {
        return !(operator == (rhs));
    }        
};


template<typename TWeight>
struct wedge_t : edge_t
{
    typedef TWeight weight_type;
    
    weight_type w;
    
    wedge_t() { }
    wedge_t(vertex_t s_, vertex_t t_, weight_type w_) : s(s_), t(t_), w(w_) { }        
};

struct no_edge_weight { };

template<typename TWeight>
struct get_edge_type
{
    typedef wedge_t<TWeight> edge_type;    
};

template<>
struct get_edge_type<no_edge_weight>
{
    typedef edge_t edge_type;
};


template<class TGraph>
inline vertex_t source(const edge_t& e, const TGraph& g)
{
    return e.s;
}

template<class TGraph>
inline vertex_t target(const edge_t& e, const TGraph& g)
{
    return e.t;
}


// basic iterators

class vertex_iterator_t 
        : public boost::iterator_facade<
                vertex_iterator_t, 
                vertex_t, 
                boost::bidirectional_traversal_tag,
                vertex_t>
{
public:
    vertex_iterator_t(graph_size_t i) : _i(i)
    {
    }
    
    vertex_t dereference() const
    {
        return _i;
    }
    
    bool equal(const vertex_iterator& r) const
    {
        return _i == r._i;
    }
    
    void increment()
    {
        ++ _i;
    }
    
    void decrement()
    {
        -- _i;
    }    
    
private:
    graph_size_t _i;
};


class vertex_array_iterator_t : 
    public boost::iterator_facade<
            vertex_array_iterator_t, 
            vertex_t, 
            boost::bidirectional_traversal_tag,
            const vertex_t&>
{
public:
    vertex_array_iterator_t(const vertex_t *pv) : _pv(pv) { }
    
    const vertex_t& dereference() const
    {
        return *_pv;
    }
    
    bool equal(const vertex_array_iterator& r) const
    {
        return _pv == r._pv;
    }
    
    void increment()
    {
        ++ _pv;
    }
    
    void decrement()
    {
        -- _pv;
    }    
    
private:
    const vertex_t *_pv;
    
};


class adj_edge_iterator_t 
        : public boost::iterator_facade<
                adj_edge_iterator_t, 
                edge_t, 
                boost::bidirectional_traversal_tag, 
                edge_t>        
{
public:
    adj_edge_iterator_t(vertex_t s, const vertex_t *ts) : _s(s), _ts(ts) { }
    
    edge_t dereference() const
    {
        return edge_t(_s, *_ts);
    }
    
    bool equal(const adj_edge_iterator& r) const
    {
        return _s == r._s && _ts == r._ts;
    }
    
    void increment() 
    {
        ++ _ts;
    }
    
    void decrement()
    {
        -- _ts;
    }
            
private:
    vertex_t _s;
    const vertex_t *_ts;        
};


// Graph classes
    
template<typename TWeight, typename TCategory=boost::incidence_graph_tag>        
class CRefAdjList
{
public:
    typedef vertex_t                                    vertex_type;
    typedef typename get_edge_type<TWeight>::edge_type  edge_type;
        
    // BGL-compliant definitions
    
    typedef vertex_type     vertex_descriptor;
    typedef edge_type       edge_descriptor;
    typedef graph_size_t    vertices_size_type;
    typedef graph_size_t    edges_size_type;
    typedef graph_size_t    degree_size_type;
    
    typedef TCategory traversal_category;
    typedef boost::directed_tag directed_category;
    typedef boost::disallow_parallel_edge_tag edge_parallel_category;           
    
    typedef vertex_iterator_t       vertex_iterator;    
    typedef adj_edge_iterator_t     out_edge_iterator;
    typedef vertex_array_iterator_t adjacency_iterator;  
    
    typedef std::pair<vertex_iterator, vertex_iterator> vertex_iterator_pair;
    typedef std::pair<out_edge_iterator, out_edge_iterator> out_edge_iterator_pair;
    typedef std::pair<adjacency_iterator, adjacency_iterator> adjacency_iterator_pair;
    
public:
    
    CRefAdjList(graph_size_t n, graph_size_t m, 
            const graph_size_t *outdegs, const graph_size_t *offsets, const vertex_t* adj_vertices)
    : m_num_vertices(n), m_num_vertices(m)
    , m_degress(outdegs), m_offsets(offsets), m_adj_vertices(adj_vertices)
    {
    }
    
    
    graph_size_t num_vertices() const
    {
        return m_num_vertices;
    }
    
    graph_size_t num_edges() const
    {
        return m_num_edges;
    }
        
    vertex_iterator vertices_begin() const
    {
        return 0;
    }
    
    vertex_iterator vertices_end() const
    {
        return m_num_vertices;
    }
    
    graph_size_t out_degree(vertex_t u) const
    {
        return m_degrees[u];
    }
    
    adjacency_iterator adj_vertices_begin(vertex_t u) const
    {
        return _adj_begin(u);
    }
    
    adjacency_iterator adj_vertices_end(vertex_t u) const
    {
        return _adj_end(u);
    }
    
    out_edge_iterator out_edges_begin(vertex_t u) const
    {
        return out_edge_iterator(u, _adj_begin(u));
    }
    
    out_edge_iterator out_edges_end(vertex_t u) const
    {
        return out_edge_iterator(u, _adj_end(u));
    }
    
private:
    const vertex_t *_adj_begin(vertex_t u) const
    {
        return m_adj_vertices + m_offsets[u];
    }
    
    const vertex_t *_adj_end(vertex_t u) const
    {
        return m_adj_vertices + (m_offsets[u] + m_degrees[u]);
    }
    
    
private:
    graph_size_t    m_num_vertices;
    graph_size_t    m_num_edges;    
    
    graph_size_t    *m_degrees;
    graph_size_t    *m_offsets;     
    vertex_t        *m_adj_vertices;
    
}; // end class AdjList



// BGL-compliant adapting functions


template<typename TWeight, typename TCategory>
inline graph_size_t num_vertices(const CRefAdjList<TWeight, TCategory>& g)
{
    return g.num_vertices();
}

template<typename TWeight, typename TCategory>
inline CRefAdjList<TWeight, TCategory>::vertex_iterator_pair
vertices(const CRefAdjList<TWeight, TCategory>& g)
{
    return std::make_pair(g.vertices_begin(), g.vertices_end());
}

template<typename TWeight, typename TCategory>
inline graph_size_t num_edges(const CRefAdjList<TWeight, TCategory>& g)
{
    return g.num_edges();
}



template<class TGraph>
inline graph_size_t out_degree(vertex_t u, const TGraph& g)
{
    return g.out_degree(u);
}

template<class TGraph>
inline graph_size_t in_degree(vertex_t u, const TGraph& g)
{
    return g.in_degree(u);
}

template<class TGraph>
inline graph_size_t degree(vertex_t u, const TGraph& g)
{
    return g.degree(u);
}

template<class TGraph>
inline TGraph::adjacency_iterator_pair 
adjacent_vertices(vertex_t u, const TGraph& g)
{
    return std::make_pair(g.adj_vertices_begin(u), g.adj_vertices_end(u));
}

template<class TGraph>
inline TGraph::out_edge_iterator_pair
out_edges(vertex_t u, const TGraph& g)
{
    return std::make_pair(g.out_edges_begin(u), g.out_edges_end(u));
}

template<class TGraph>
inline TGraph::in_edge_iterator_pair
in_edges(vertex_t u, const TGraph& g)
{
    return std::make_pair(g.in_edges_begin(u), g.in_edges_end(u));
}





}


#endif
