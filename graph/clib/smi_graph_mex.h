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
#include "smi_graph.h"

namespace smi
{    
    
    
inline DGraphSpec get_dgraph_spec(bcs::matlab::const_marray mGr)
{
    DGraphSpec s;
    
    char dty = (char)mGr.get_field(0, "dty").get_scalar<mxChar>();
    s.n = mGr.get_field(0, "n").get_scalar<gindex_t>();
    s.m = mGr.get_field(0, "m").get_scalar<gindex_t>();
    s.edges = mGr.get_field(0, "edges").data<DEdge>();
      
    if (dty == 'u') s.m *= 2;
                
    bool has_nbs = mGr.get_field(0, "has_nbs").get_scalar<bool>();
    
    if (has_nbs)
    {
        s.out_nbs = mGr.get_field(0, "o_nbs").data<gindex_t>();
        s.out_edges = mGr.get_field(0, "o_eds").data<gindex_t>();
        s.out_degs = mGr.get_field(0, "o_degs").data<gindex_t>();
        s.offsets = mGr.get_field(0, "o_os").data<gindex_t>();
    }
    else
    {
        s.out_nbs = 0;
        s.out_edges = 0;
        s.out_degs = 0;
        s.offsets = 0;
    }
    
    return s;
}
        
    
}
