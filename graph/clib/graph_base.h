/********************************************************************
 *
 *  graph_base.h
 *
 *  The base facilities for graph library
 *  
 *  Created by Dahua Lin, on Oct 8, 2010
 *
 ********************************************************************/

#ifndef SMI_GRAPH_BASE_H
#define SMI_GRAPH_BASE_H

#include <utility>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_concepts.hpp>

#include <valarray>


namespace smi
{

/************************************************
 *
 *  Basic types
 *
 ************************************************/
    
typedef int graph_size_t;        

struct vertex_t
{
    graph_size_t i;
    
    vertex_t( ) { }
    vertex_t(graph_size_t i_) : i(i_) { }
    
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
    graph_size_t i;
    
    edge_t( ) { }
    edge_t(graph_size_t i_) : i(i_) { }
    
    bool operator == (const edge_t& rhs) const
    {
        return i == rhs.i;
    }
    
    bool operator != (const edge_t& rhs) const
    {
        return i != rhs.i;
    }
};


struct no_edge_weight { };


struct edge_info
{
    vertex_t s;
    vertex_t t;
    
    edge_info( ) { }
    edge_info(vertex_t s_, vertex_t t_) : s(s_), t(t_) { }
};


template<typename TWeight>
struct wedge_info : public edge_info
{
    typedef TWeight weight_type;
    
    weight_type w;
    
    wedge_info() { }
    wedge_info(vertex_t s_, vertex_t t_, weight_type w_) 
    : edge_info(s_, t_), w(w_) { }
    
    template<typename TGraph>
    static wedge_info<TWeight> from_edge(const edge_t& e, const TGraph& g)
    {
        return wedge_info<TWeight>(g.get_source(e), g.get_target(e), g.get_weight(e));
    }
};


template<>
struct wedge_info<no_edge_weight> : public edge_info
{
    typedef no_edge_weight weight_type;
    
    wedge_info() { }
    wedge_info(vertex_t s_, vertex_t t_) 
    : edge_info(s_, t_) { }
    
    template<typename TGraph>
    static wedge_info<no_edge_weight> from_edge(const edge_t& e, const TGraph& g)
    {
        return wedge_info<no_edge_weight>(g.get_source(e), g.get_target(e));
    }
};



/************************************************
 *
 *  Basic iterators
 *
 ************************************************/

class vertex_iterator_t 
        : public boost::iterator_facade<
                vertex_iterator_t, 
                vertex_t, 
                boost::bidirectional_traversal_tag,
                vertex_t>
{
public:
    vertex_iterator_t() { }
    
    vertex_iterator_t(graph_size_t i) : _i(i)
    {
    }
    
    vertex_t dereference() const
    {
        return _i;
    }
    
    bool equal(const vertex_iterator_t& r) const
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



class edge_iterator_t 
        : public boost::iterator_facade<
                edge_iterator_t, 
                edge_t, 
                boost::bidirectional_traversal_tag,
                edge_t>
{
public:
    edge_iterator_t() { }
    
    edge_iterator_t(graph_size_t i) : _i(i)
    {
    }
    
    edge_t dereference() const
    {
        return _i;
    }
    
    bool equal(const edge_iterator_t& r) const
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
    vertex_array_iterator_t() { }
    
    vertex_array_iterator_t(const vertex_t *p) : _p(p) { }
    
    const vertex_t& dereference() const
    {
        return *_p;
    }
    
    bool equal(const vertex_array_iterator_t& r) const
    {
        return _p == r._p;
    }
    
    void increment()
    {
        ++ _p;
    }
    
    void decrement()
    {
        -- _p;
    }    
    
private:
    const vertex_t *_p;
    
};



class edge_array_iterator_t : 
    public boost::iterator_facade<
            edge_array_iterator_t, 
            edge_t, 
            boost::bidirectional_traversal_tag,
            const edge_t&>
{
public:
    edge_array_iterator_t() { }
    
    edge_array_iterator_t(const edge_t *p) : _p(p) { }
    
    const edge_t& dereference() const
    {
        return *_p;
    }
    
    bool equal(const edge_array_iterator_t& r) const
    {
        return _p == r._p;
    }
    
    void increment()
    {
        ++ _p;
    }
    
    void decrement()
    {
        -- _p;
    }    
    
private:
    const edge_t *_p;
    
};



/************************************************
 *
 *  Property map classes
 *
 ************************************************/

// for vertex index

class vertex_index_dmap
{
public:
    typedef vertex_t key_type;
    typedef graph_size_t value_type;
    typedef graph_size_t reference;
    typedef boost::readable_property_map_tag category;
    
public:
    reference get(const key_type& v) const
    {
        return v.i;
    }    
};

inline graph_size_t get(const vertex_index_dmap& m, const vertex_t& v)
{
    return m.get(v);
}

template<typename TGraph>
inline vertex_index_dmap get(boost::vertex_index_t, const TGraph& g)
{
    return vertex_index_dmap();
}

template<typename TGraph>
inline graph_size_t get(boost::vertex_index_t, const TGraph& g, const vertex_t &v)
{
    return v.i;
}


// for edge index

class edge_index_dmap
{
public:
    typedef edge_t key_type;
    typedef graph_size_t value_type;
    typedef graph_size_t reference;
    typedef boost::readable_property_map_tag category;
    
public:
    reference get(const key_type& e) const
    {
        return e.i;
    }    
};

inline graph_size_t get(const edge_index_dmap& m, const edge_t& e)
{
    return m.get(e);
}

template<typename TGraph>
inline edge_index_dmap get(boost::edge_index_t, const TGraph& g)
{
    return edge_index_dmap();
}

template<typename TGraph>
inline graph_size_t get(boost::edge_index_t, const TGraph& g, const edge_t &e)
{
    return e.i;
}


// for edge weight

template<typename TWeight>
class edge_weight_crefmap
{
public:
    typedef edge_t key_type;
    typedef TWeight value_type;
    typedef const TWeight& reference;    
    typedef boost::readable_property_map_tag category;
            
public:
    edge_weight_crefmap(const TWeight *ws) : m_ws(ws) { }
    
    reference get(const key_type& e) const
    {
        return m_ws[e.i];
    }
    
private:
    const TWeight *m_ws;    
};

template<typename TWeight>
inline const TWeight& get(const edge_weight_crefmap<TWeight>& wmap, const edge_t& e)
{
    return wmap.get(e);
}



/************************************************
 *
 *  Functions to extract attributes 
 *  of vertices / edges from a graph
 *
 ************************************************/

template<class TGraph>
inline vertex_t source(const edge_t& e, const TGraph& g)
{
    return g.get_source(e.i);
}

template<class TGraph>
inline vertex_t target(const edge_t& e, const TGraph& g)
{
    return g.get_target(e.i);
}

template<class TGraph>
inline graph_size_t out_degree(const vertex_t& u, const TGraph& g)
{
    return g.out_degree(u);
}

template<class TGraph>
inline graph_size_t in_degree(const vertex_t& v, const TGraph& g)
{
    return g.in_degree(v);
}

template<class TGraph>
inline graph_size_t degree(const vertex_t& v, const TGraph& g)
{
    return g.degree(v);
}


template<class TGraph>
inline std::pair<typename TGraph::out_edge_iterator, typename TGraph::out_edge_iterator> 
out_edges(const vertex_t& u, const TGraph& g)
{
    return std::make_pair(g.out_edges_begin(u), g.out_edges_end(u));
}

template<class TGraph>
inline std::pair<typename TGraph::in_edge_iterator, typename TGraph::in_edge_iterator> 
in_edges(const vertex_t& v, const TGraph& g)
{
    return std::make_pair(g.in_edges_begin(v), g.in_edges_end(v));
}


template<class TGraph>
inline std::pair<typename TGraph::adjacency_iterator, typename TGraph::adjacency_iterator> 
adjacent_vertices(const vertex_t& u, const TGraph& g)
{
    return std::make_pair(g.adj_vertices_begin(u), g.adj_vertices_end(u));
}


/************************************************
 *
 *  Property maps (based on valarray)
 *
 ***********************************************/

template<typename T>
class VertexRefMap
{
public:
    typedef vertex_t key_type;
    typedef T value_type;
    typedef value_type& reference;    
    typedef boost::lvalue_property_map_tag category;    
    
    VertexRefMap(T *vals) : m_values(vals)
    {
    }
                 
    const value_type& operator[] (const vertex_t& v) const
    {
        return m_values[v.i];
    }
    
    value_type& operator[] (const vertex_t& v)
    {
        return m_values[v.i];
    }
    
private:
    T *m_values;    
};


template<typename T>
const T& get(const VertexRefMap<T>& m, const vertex_t& v)
{
    return m[v.i];
}

template<typename T>
const T& get(const VertexRefMap<T>& m, unsigned long i)
{
    return m[i];
}

template<typename T>
void put(VertexRefMap<T>& m, const vertex_t& v, const T& c)
{
    m[v.i] = c;
}

template<typename T>
void put(VertexRefMap<T>& m, unsigned long i, const T& c)
{
    m[i] = c;
}


/************************************************
 *
 *  A simple color map for use in some algorithms
 *
 ************************************************/

struct GraphSearchColorValue
{
    char _c;
    
    GraphSearchColorValue() { }
    
    GraphSearchColorValue(char c) : _c(c) { }
               
    bool operator == (const GraphSearchColorValue& r) const
    {
        return _c == r._c;
    }
    
    bool operator != (const GraphSearchColorValue& r) const
    {
        return _c != r._c;
    }    
    
    bool is_white() const { return _c == char(0); }
    bool is_gray()  const { return _c == char(1); }
    bool is_black() const { return _c == char(2); }
    
    static GraphSearchColorValue white() { return char(0); }
    static GraphSearchColorValue gray()  { return char(1); }
    static GraphSearchColorValue black() { return char(2); }
};


}




#endif


