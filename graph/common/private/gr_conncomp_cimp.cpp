/********************************************************************
 *
 *  gr_conncomp_cimp.cpp
 *
 *  The C++ mex implementation of gr_conncomp
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 *******************************************************************/


#include "../../clib/graph_mex.h"

#include <vector>
#include <valarray>

#include <boost/graph/connected_components.hpp>

typedef boost::default_color_type color_t;

using namespace smi;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    matlab_graph_repr gr(prhs[0]);    
    CRefAdjList<no_edge_weight, boost::undirected_tag> g = gr.to_cref_adjlist_ud();
    
    graph_size_t n = num_vertices(g);
        
    std::valarray<color_t> colors(boost::white_color, n);
    std::valarray<graph_size_t> labels(n);
        
    VertexRefMap<color_t> cmap = &(colors[0]);
    VertexRefMap<graph_size_t> lmap = &(labels[0]);
    
    boost::connected_components(g, lmap, boost::color_map(cmap));
    
    plhs[0] = iter_to_matlab_row(&(labels[0]), labels.size(), vertex_to_mindex());
}

