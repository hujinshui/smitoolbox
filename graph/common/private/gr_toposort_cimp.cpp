/********************************************************************
 *
 *  gr_toposort_cimp.cpp
 *
 *  The C++ mex implementation of topological sorting
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 *******************************************************************/


#include "../../clib/graph_mex.h"

#include <vector>
#include <valarray>
#include <algorithm>
#include <iterator>

#include <boost/graph/topological_sort.hpp>

typedef boost::default_color_type color_t;

using namespace smi;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    matlab_graph_repr gr(prhs[0]);    
    CRefAdjList<no_edge_weight> g = gr.to_cref_adjlist();
    
    graph_size_t n = num_vertices(g);
    
    std::vector<vertex_t> tord;
    tord.reserve(n);    
    std::valarray<color_t> colors(n);
    
    VertexRefMap<color_t> cmap = &(colors[0]);
    
    bool is_dag = true;            
    try
    {
        using boost::color_map;
        
        boost::topological_sort(g, std::back_inserter(tord), color_map(cmap));
    }
    catch( boost::not_a_dag )
    {
        is_dag = false;
    }
    
    if (is_dag)
    {
        std::reverse(tord.begin(), tord.end());
        plhs[0] = iter_to_matlab_row(tord.begin(), tord.size(), vertex_to_mindex());        
    }
    else
    {
        plhs[0] = create_matlab_scalar<int>(-1);
    }    
    
}

