/********************************************************************
 *
 *  smi_graph_mex.h
 *
 *  The mex interface for SMI graphs
 *
 *  Created by Dahua Lin, on Oct 28, 2011
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/graph/gedgelist_view.h>
#include <bcslib/graph/ginclist_view.h>

namespace smi
{

// Typedefs    
    
typedef bcs::int32_t gint;
typedef bcs::gvertex<gint> vertex_t;
typedef bcs::gedge<gint> edge_t;
typedef bcs::gvertex_pair<gint> vertex_pair;

typedef bcs::gedgelist_view<gint> gedgelist_t;
typedef bcs::ginclist_view<gint> ginclist_t;

typedef bcs::natural_vertex_iterator<gint>::type natural_viter_t;
typedef bcs::natural_edge_iterator<gint>::type natural_eiter_t;

// Interoperability

bool isdirected_from_char(char dty)
{
    return dty == 'd' || dty == 'D';
}


inline gedgelist_t mx_to_gedgelist(bcs::matlab::const_marray mG)
{
    char dty = (char)mG.get_field(0, "dty").get_scalar<mxChar>();
    gint n = mG.get_field(0, "n").get_scalar<gint>();
    gint m = mG.get_field(0, "m").get_scalar<gint>();
    const vertex_pair* edges = mG.get_field(0, "edges").data<vertex_pair>();
    
    bool is_dir = isdirected_from_char(dty);
    
    return gedgelist_t(n, m, is_dir, edges);
}

inline ginclist_t mx_to_ginclist(bcs::matlab::const_marray mG)
{    
    char dty = (char)mG.get_field(0, "dty").get_scalar<mxChar>();
    gint n = mG.get_field(0, "n").get_scalar<gint>();
    gint m = mG.get_field(0, "m").get_scalar<gint>();
    const vertex_pair* edges = mG.get_field(0, "edges").data<vertex_pair>(); 
    
    const vertex_t* nbs = mG.get_field(0, "o_nbs").data<vertex_t>();
    const edge_t* eds = mG.get_field(0, "o_eds").data<edge_t>();
    const gint* degs = mG.get_field(0, "o_degs").data<gint>();
    const gint* ofs = mG.get_field(0, "o_os").data<gint>();
    
    bool is_dir = isdirected_from_char(dty);
    
    return ginclist_t(n, m, is_dir, edges, nbs, eds, degs, ofs);
}


} // end namespace smi

