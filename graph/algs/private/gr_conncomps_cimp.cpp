/***************************************************************
 *
 *  gr_conncomps_cimp.cpp
 *
 *  The C++ mex implementation of gr_conncomps_cimp
 *
 *  Created by Dahua Lin, on Jan 26, 2012
 *
 ***************************************************************/

#include "../../clib/smi_graph_mex.h"
#include <bcslib/graph/graph_traversal.h>

#include <vector>
#include <list>

using namespace bcs;
using namespace bcs::matlab;
using namespace smi;

struct ccs_recorder
{
    ccs_recorder() : cur_comp(0) { }
    
    void new_component() 
    {
        comp_list.push_back(std::vector<vertex_t>());
        cur_comp = &(comp_list.back());
    }
    
    void end_component()
    {
    }
    
    void add_vertex(const vertex_t& v)
    {
        cur_comp->push_back(v);
    }
    
    std::list<std::vector<vertex_t> > comp_list;
    std::vector<vertex_t>* cur_comp;
    
    marray moutput() const
    {
        index_t nc = (index_t)comp_list.size();
        marray mCCS = create_mcell_array(1, nc);
        
        std::list<std::vector<vertex_t> >::const_iterator p =
                comp_list.begin();
        
        for (index_t i = 0; i < nc; ++i)
        {
            const std::vector<vertex_t>& vec = *(p++);
            
            caview1d<gint> a((const gint*)(&(vec[0])), (index_t)vec.size());
            mCCS.set_cell(i, to_matlab_row(a));
        }
        
        return mCCS;
    }
};


/**
 * The main entry
 *
 * Input
 *   [0] G:         the graph struct
 *  
 * Output
 *   [0] ccs:       the cell array of components
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mG(prhs[0]);
    ginclist_t G = mx_to_ginclist(mG);
    
    // solve
    
    ccs_recorder rec;
    find_connected_components(G, rec);
    
    // output
    
    marray mCCS = rec.moutput();
    plhs[0] = mCCS.mx_ptr();
}

BCSMEX_MAINDEF
