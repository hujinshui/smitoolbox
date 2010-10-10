/********************************************************************
 *
 *  graph_adjlist.h
 *
 *  The BGL-compliant adjacency list classes which are designed 
 *  for efficient interaction with MATLAB representation.
 *  
 *  Created by Dahua Lin, on Oct 8, 2010
 *
 ********************************************************************/


#ifndef SMI_GRAPH_ADJLIST_H
#define SMI_GRAPH_ADJLIST_H


#include "graph_base.h"


namespace smi
{

/************************************************
 *
 *  Edge List classes
 *
 ************************************************/
    
    
/**
 * The graph class based on referenced edge list
 *
 * It implements the following concepts in BGL:
 *
 *      VertexListGraph
 *      EdgeListGraph
 *      
 *      PropertyGraph (with readable weightmap)
 *
 */
template<typename TWeight=no_edge_weight, 
        typename TCategory=boost::edge_list_graph_tag,
        typename TDirected=boost::directed_tag>
class CRefEdgeList
{
public:
    typedef vertex_t    vertex_type;
    typedef edge_t      edge_type;
    typedef TWeight     edge_weight_type;
    
    typedef wedge_info<TWeight> edge_info_type;
        
    typedef vertex_iterator_t   vertex_iterator;
    typedef edge_iterator_t     edge_iterator; 
            
public:
    
    CRefEdgeList(graph_size_t n, graph_size_t m, 
            const vertex_type *srcs, 
            const vertex_type *tars, 
            const edge_weight_type *ws)
    : m_num_vertices(n), m_num_edges(m)
    , m_sources(srcs), m_targets(tars), m_weights(ws)
    {
        
    }
    
public:
    
    graph_size_t num_vertices() const
    {
        return m_num_vertices;
    }
    
    graph_size_t num_edges() const
    {
        return m_num_edges;
    }
    
    const vertex_type *sources() const
    {
        return m_sources;
    }
    
    const vertex_type *targets() const
    {
        return m_targets;
    }
    
    const edge_weight_type *weights() const
    {
        return m_weights;
    }
    
    const vertex_type& get_source(const edge_type& e) const
    {
        return m_sources[e.i];
    }
    
    const vertex_type& get_target(const edge_type& e) const
    {
        return m_targets[e.i];
    }
    
    const edge_weight_type& get_weight(const edge_type& e) const
    {
        return m_weights[e.i];
    }
    
public:
    
    vertex_iterator vertices_begin() const
    {
        return 0;
    }
    
    vertex_iterator vertices_end() const
    {
        return m_num_vertices;
    }
    
    edge_iterator edges_begin() const
    {
        return 0;
    }
    
    edge_iterator edges_end() const
    {
        return m_num_edges;
    }           
    
private:
    graph_size_t m_num_vertices;
    graph_size_t m_num_edges;
    
    const vertex_type *m_sources;
    const vertex_type *m_targets;
    const edge_weight_type *m_weights;
    
}; // end class CRefEdgeList


// BGL-compliant adapting functions for CRefEdgeList

template<typename TWeight, typename TCategory, typename TDirected>
inline graph_size_t
num_vertices(const CRefEdgeList<TWeight, TCategory, TDirected>& g)
{
    return g.num_vertices();
}


template<typename TWeight, typename TCategory, typename TDirected>
inline std::pair<vertex_iterator_t, vertex_iterator_t> 
vertices(const CRefEdgeList<TWeight, TCategory, TDirected>& g)
{
    return std::make_pair(g.vertices_begin(), g.vertices_end());
}


template<typename TWeight, typename TCategory, typename TDirected>
inline graph_size_t
num_edges(const CRefEdgeList<TWeight, TCategory, TDirected>& g)
{
    return g.num_edges();
}


template<typename TWeight, typename TCategory, typename TDirected>
inline std::pair<edge_iterator_t, edge_iterator_t> 
edges(const CRefEdgeList<TWeight, TCategory, TDirected>& g)
{
    return std::make_pair(g.edges_begin(), g.edges_end());
}


template<typename TWeight, typename TCategory, typename TDirected>
inline vertex_index_dmap
get(boost::vertex_index_t, const CRefEdgeList<TWeight, TCategory, TDirected>& g)
{
    return vertex_index_dmap();
}

template<typename TWeight, typename TCategory, typename TDirected>
inline vertex_index_dmap
get(boost::edge_index_t, const CRefEdgeList<TWeight, TCategory, TDirected>& g)
{
    return edge_index_dmap();
}

template<typename TWeight, typename TCategory, typename TDirected>
inline edge_weight_crefmap<TWeight>
get(boost::edge_weight_t, const CRefEdgeList<TWeight, TCategory, TDirected>& g)
{
    return g.weights();
}



/**
 *  Adjacency List classes
 *
 *  As a derived class of CRefEdgeList, it 
 *  implements
 *
 *      VertexListGraph
 *      EdgeListGraph
 *      PropertyGraph (with readable weightmap)
 *
 *  In addition, it implements the following 
 *  concept:
 *
 *      IncidenceGraph
 *      AdjacencyGraph
 */
template<typename TWeight=no_edge_weight, 
        typename TCategory=boost::incidence_graph_tag,
        typename TDirected=boost::directed_tag>
class CRefAdjList : public CRefEdgeList<TWeight, TCategory, TDirected>
{
public:
    typedef CRefEdgeList<TWeight, TCategory, TDirected> _base_class;
    
    typedef typename _base_class::vertex_type       vertex_type;
    typedef typename _base_class::edge_type         edge_type;
    typedef typename _base_class::edge_weight_type  edge_weight_type;
            
    typedef typename _base_class::vertex_iterator   vertex_iterator;
    typedef typename _base_class::edge_iterator     edge_iterator;
    
    typedef edge_iterator_t         out_edge_iterator;
    typedef vertex_array_iterator_t adjacency_iterator;
       
public:
    CRefAdjList(graph_size_t n, graph_size_t m, 
            const vertex_type *srcs, 
            const vertex_type *tars, 
            const edge_weight_type *ws,
            const graph_size_t *out_degs, 
            const graph_size_t *out_offsets)
    : _base_class(n, m, srcs, tars, ws)
    , m_out_degs(out_degs), m_out_offsets(out_offsets)
    {        
    }
    
    
public:
    
    graph_size_t out_degree(const vertex_type& u) const
    {
        return m_out_degs[u.i];
    }
    
    out_edge_iterator out_edges_begin(const vertex_type& u) const
    {
        return adj_begin(u);
    }
    
    out_edge_iterator out_edges_end(const vertex_type& u) const
    {
        return adj_end(u);
    }
    
    adjacency_iterator adj_vertices_begin(const vertex_type& u) const
    {
        return this->targets() + adj_begin(u);
    }
    
    adjacency_iterator adj_vertices_end(const vertex_type& u) const
    {
        return this->targets() + adj_end(u);
    }
    
private:
    graph_size_t adj_begin(const vertex_type& u) const
    {
        return m_out_offsets[u.i];
    }
    
    graph_size_t adj_end(const vertex_type& u) const
    {
        return m_out_offsets[u.i] + m_out_degs[u.i];
    }
    
    
private:
    const graph_size_t *m_out_degs;
    const graph_size_t *m_out_offsets;
        
    
};  // end class CRefAdjList

}


namespace boost
{

template<typename TWeight, typename TCategory, typename TDirected>
struct graph_traits<smi::CRefEdgeList<TWeight, TCategory, TDirected> >
{
    typedef smi::CRefEdgeList<TWeight, TCategory, TDirected> graph_type;
    
    typedef typename graph_type::vertex_type    vertex_descriptor;
    typedef typename graph_type::edge_type      edge_descriptor;
    
    typedef smi::graph_size_t   vertices_size_type;
    typedef smi::graph_size_t   edges_size_type;
    
    typedef TCategory traversal_category;
    typedef TDirected directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category;
    
    typedef typename graph_type::vertex_iterator vertex_iterator;
    typedef typename graph_type::edge_iterator edge_iterator;
};


template<typename TWeight, typename TCategory, typename TDirected>
struct property_map<smi::CRefEdgeList<TWeight, TCategory, TDirected>, boost::vertex_index_t>
{    
    typedef smi::vertex_index_dmap type;
    typedef smi::vertex_index_dmap const_type;
};


template<typename TWeight, typename TCategory, typename TDirected>
struct graph_traits<smi::CRefAdjList<TWeight, TCategory, TDirected> >
{
    typedef smi::CRefAdjList<TWeight, TCategory, TDirected> graph_type;
    
    typedef typename graph_type::vertex_type    vertex_descriptor;
    typedef typename graph_type::edge_type      edge_descriptor;
    
    typedef smi::graph_size_t   vertices_size_type;
    typedef smi::graph_size_t   edges_size_type;
    typedef smi::graph_size_t   degree_size_type;
    
    typedef TCategory traversal_category;
    typedef TDirected directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category;
    
    typedef typename graph_type::vertex_iterator vertex_iterator;
    typedef typename graph_type::edge_iterator edge_iterator;
    typedef typename graph_type::out_edge_iterator out_edge_iterator;
    typedef typename graph_type::adjacency_iterator adjacency_iterator;
};


}




#endif




