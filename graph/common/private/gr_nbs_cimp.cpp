/********************************************************************
 *
 *  gr_nbs_cimp.cpp
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
    gedgelist_t E = mx_to_gedgelist(mG);
    
    // prepare storage
    
    gint n = E.nvertices();
    gint ma = E.is_directed() ? E.nedges() : 2 * E.nedges();
    
    marray mNbs = create_marray<gint>(1, ma);
    marray mEds = create_marray<gint>(1, ma);
    marray mDegs = create_marray<gint>(1, n);
    marray mOfs = create_marray<gint>(1, n);
    
    vertex_t *nbs = mNbs.data<vertex_t>();
    edge_t *eds = mEds.data<edge_t>();
    gint *degs = mDegs.data<gint>();
    gint *ofs = mOfs.data<gint>();
    
    // build neighborhood
    
    gint *temp = new gint[n];
    
    const vertex_pair* vps = E.vertex_pairs_begin();
    prepare_ginclist_arrays(n, ma, vps, nbs, eds, degs, ofs, temp);
    
    delete[] temp;
    
    // export
        
    plhs[0] = mNbs.mx_ptr();
    plhs[1] = mEds.mx_ptr();
    plhs[2] = mDegs.mx_ptr();
    plhs[3] = mOfs.mx_ptr();        
}


BCSMEX_MAINDEF

