/********************************************************************
 *
 *  gr_conncomp_cimp.cpp
 *
 *  The C++ mex implementation of gr_conncomp
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 *******************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/matlab/mgraph.h>
#include <bcslib/graph/bgl_port.h>

#include <boost/graph/connected_components.hpp>

#include <vector>
#include <valarray>

using namespace bcs;
using namespace bcs::matlab;


void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_mgraph mG(prhs[0]);    
    gr_adjlist<gr_undirected> g = to_gr_adjlist<gr_undirected>(mG);
    
    gr_size_t n = num_vertices(g);
    
    // main
    
    typedef boost::default_color_type color_t;
    std::valarray<color_t> colors(boost::white_color, n);
    std::valarray<gr_size_t> labels(n);
        
    vertex_ref_map<color_t> cmap = &(colors[0]);
    vertex_ref_map<gr_size_t> lmap = &(labels[0]);
    
    boost::connected_components(g, lmap, boost::color_map(cmap));
    
    // output
    
    marray mL = create_marray<int32_t>(1, n);
    int32_t *L = mL.data<int32_t>();
    for (gr_size_t i = 0; i < n; ++i)
    {
        L[i] = (int32_t)(labels[i]) + 1;
    }
    
    plhs[0] = mL.mx_ptr();
   
}


BCSMEX_MAINDEF

