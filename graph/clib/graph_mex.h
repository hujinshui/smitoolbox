/********************************************************************
 *
 *  graph_mex.h
 *
 *  The mex interface for graph classes
 *
 *  Created by Dahua Lin, on Oct 8, 2010
 *
 ********************************************************************/


#ifndef SMI_GRAPH_MEX_H
#define SMI_GRAPH_MEX_H


#include "../../base/clib/matlab_types.h"
#include "../../base/clib/marray.h"
#include "../../base/clib/functors.h"

#include "graph_adjlist.h"


namespace smi
{

struct matlab_graph_repr
{
    MArray mG;
    
    matlab_graph_repr(const MArray& G) : mG(G) { }
    
    // field extraction
    
    bool has_weight() const
    {
        return !mG.get_property("ew").is_empty();
    }
    
    mxClassID weight_class() const
    {
        return mG.get_property("ew").class_id();
    }
    
    graph_size_t n() const
    {
        return (graph_size_t)mG.get_property("nv").get_double_scalar();
    }
    
    graph_size_t m() const
    {
        return (graph_size_t)mG.get_property("ne").get_double_scalar();
    }
    
    const vertex_t* s() const
    {
        return mG.get_property("es").get_data<vertex_t>();
    }
    
    const vertex_t *t() const
    {
        return mG.get_property("et").get_data<vertex_t>();
    }
    
    template<typename TWeight>
    const TWeight *w() const
    {
        MArray mx = mG.get_property("ew");
        return mx.get_data<TWeight>();
    }
    
    const graph_size_t *o_degs() const
    {
        return mG.get_property("o_ds").get_data<graph_size_t>();
    }
    
    const graph_size_t *o_offsets() const
    {
        return mG.get_property("o_os").get_data<graph_size_t>();
    }
    
    const edge_t *o_edges() const
    {
        return mG.get_property("o_es").get_data<edge_t>();
    }
    
    const vertex_t *o_neighbors() const
    {
        return mG.get_property("o_ns").get_data<vertex_t>();
    }
    

    // conversion to ref graph class
    
    inline CRefEdgeList<no_edge_weight> to_cref_edgelist() const
    {    
        return CRefEdgeList<no_edge_weight>(n(), m(), s(), t(), 0);
    } 
    
    template<typename TWeight>
    inline CRefEdgeList<TWeight> to_cref_wedgelist() const
    {
        return CRefEdgeList<TWeight>(
                n(), m(), s(), t(), w<TWeight>());
    } 
    
    
    template<typename TWeight>
    inline CRefEdgeList<TWeight, boost::undirected_tag> to_cref_wedgelist_ud() const
    {
        return CRefEdgeList<TWeight, boost::undirected_tag>(
                n(), m(), s(), t(), w<TWeight>());
    }    
    
    
    inline CRefAdjList<no_edge_weight> to_cref_adjlist() const
    {
        return CRefAdjList<no_edge_weight>(
                n(), m(), s(), t(), 0, 
                o_degs(), o_offsets(), o_edges(), o_neighbors());
    }    
    
    
    inline CRefAdjList<no_edge_weight, boost::undirected_tag> to_cref_adjlist_ud() const
    {
        return CRefAdjList<no_edge_weight, boost::undirected_tag>(
                n(), m(), s(), t(), 0, 
                o_degs(), o_offsets(), o_edges(), o_neighbors());
    }   
    
    template<typename TWeight>
    inline CRefAdjList<TWeight> to_cref_wadjlist() const
    {
        return CRefAdjList<TWeight>(
                n(), m(), s(), t(), w<TWeight>(), 
                o_degs(), o_offsets(), o_edges(), o_neighbors());
    }
    
    template<typename TWeight>
    inline CRefAdjList<TWeight, boost::undirected_tag> to_cref_wadjlist_ud() const
    {
        return CRefAdjList<TWeight, boost::undirected_tag>(
                n(), m(), s(), t(), w<TWeight>(), 
                o_degs(), o_offsets(), o_edges(), o_neighbors());
    }
};


// Convenient functors

struct vertex_to_index
{
    typedef const vertex_t& argument_type;
    typedef graph_size_t result_type;
    
    result_type operator() (argument_type v)
    {
        return v.i;
    }
};


struct vertex_to_mindex
{
    typedef const vertex_t& argument_type;
    typedef graph_size_t result_type;
    
    result_type operator() (argument_type v)
    {
        return v.i + 1;
    }
};


struct edge_to_index
{
    typedef const edge_t& argument_type;
    typedef graph_size_t result_type;
    
    result_type operator() (argument_type e)
    {
        return e.i;
    }
};


struct edge_to_mindex
{
    typedef const edge_t& argument_type;
    typedef graph_size_t result_type;
    
    result_type operator() (argument_type e)
    {
        return e.i + 1;
    }
};


template<typename TGraph>
struct edge_to_source_t
{
    typedef const edge_t& argument_type;
    typedef vertex_t result_type;
    
    const TGraph &_g;
    edge_to_source_t(const TGraph& g) : _g(g) { } 
    
    result_type operator() (argument_type e)
    {
        return source(e, _g);
    }    
};


template<typename TGraph>
struct edge_to_target_t
{
    typedef const edge_t& argument_type;
    typedef vertex_t result_type;
    
    const TGraph &_g;
    edge_to_target_t(const TGraph& g) : _g(g) { } 
    
    result_type operator() (argument_type e)
    {
        return target(e, _g);
    }    
};


template<typename TGraph>
struct edge_to_weight_t
{
    typedef const edge_t& argument_type;
    typedef typename TGraph::edge_weight_type result_type;
    
    const TGraph &_g;
    edge_to_weight_t(const TGraph& g) : _g(g) { } 
    
    result_type operator() (argument_type e)
    {
        return _g.get_weight(e);
    }    
};



template<typename TGraph>
inline edge_to_source_t<TGraph> edge_to_source(const TGraph& g)
{
    return edge_to_source_t<TGraph>(g);
}

template<typename TGraph>
inline edge_to_target_t<TGraph> edge_to_target(const TGraph& g)
{
    return edge_to_target_t<TGraph>(g);
}

template<typename TGraph>
inline edge_to_weight_t<TGraph> edge_to_weight(const TGraph& g)
{
    return edge_to_weight_t<TGraph>(g);
}




}


#endif

