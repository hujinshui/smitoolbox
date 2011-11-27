/********************************************************************
 *
 *  make_gr_cimp.cpp
 *
 *  The C++ mex implementation of part of make_gr functionality
 *
 *  Created by Dahua Lin, on Oct 28, 2011
 *
 ********************************************************************/


#include "../../clib/smi_graph_mex.h"

using namespace bcs;
using namespace bcs::matlab;
using namespace smi;

/**
 * The main entry
 *
 * Input
 *   [0] s:         the graph struct (without nbs)
 *  
 * Output
 *   [0] o_nbs
 *   [1] o_eds
 *   [2] o_degs
 *   [3] o_os
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const_marray mG(prhs[0]);
    DGraphSpec s = get_dgraph_spec(mG);
    
    // build the graph with neighborhood system
    
    DGraph g(s.n, s.m);
    
    for (gindex_t i = 0; i < s.m; ++i)
    {
        g.add_edge(s.edges[i]);
    }
    
    g.build_neighborhood_system();
    
    // export
    
    s = g.get_spec();
        
    plhs[0] = to_matlab_row(make_caview1d(s.out_nbs, s.m)).mx_ptr();
    plhs[1] = to_matlab_row(make_caview1d(s.out_edges, s.m)).mx_ptr();
    plhs[2] = to_matlab_row(make_caview1d(s.out_degs, s.n)).mx_ptr();
    plhs[3] = to_matlab_row(make_caview1d(s.offsets, s.n)).mx_ptr();        
}


BCSMEX_MAINDEF

