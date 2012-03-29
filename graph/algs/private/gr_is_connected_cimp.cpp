/***************************************************************
 *
 *  gr_is_connected_cimp.cpp
 *
 *  The C++ mex implementation of gr_is_connected_cimp
 *
 *  Created by Dahua Lin, on Jan 26, 2012
 *
 ***************************************************************/

#include "../../clib/smi_graph_mex.h"
#include <bcslib/graph/graph_traversal.h>

using namespace bcs;
using namespace bcs::matlab;
using namespace smi;

/**
 * The main entry
 *
 * Input
 *   [0] G:         the graph struct
 *  
 * Output
 *   [0] tf:        whether G is connected (logical scalar)
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mG(prhs[0]);
    ginclist_t G = mx_to_ginclist(mG);
    
    // solve
    
    size_t nv = (size_t)G.nvertices();
    bool tf = true;
    
    if (nv > 1)
    {
        size_t cr = count_reachable_vertices(G, make_gvertex(gint(1)));
        tf = (cr == nv - 1);
    }
    
    // output
    
    plhs[0] = mxCreateLogicalScalar((mxLogical)tf);
}

BCSMEX_MAINDEF
