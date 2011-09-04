/********************************************************************
 *
 *  gr_toposort_cimp.cpp
 *
 *  The C++ mex implementation of topological sorting
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 *******************************************************************/


#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/matlab/mgraph.h>
#include <bcslib/graph/bgl_port.h>

#include <boost/graph/topological_sort.hpp>

#include <vector>
#include <valarray>
#include <iterator>

using namespace bcs;
using namespace bcs::matlab;


void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_mgraph mG(prhs[0]);   
    gr_adjlist<gr_undirected> g = to_gr_adjlist<gr_undirected>(mG);
    
    gr_size_t n = num_vertices(g);
    
    // prepare
    
    typedef boost::default_color_type color_t;
    std::vector<vertex_t> tord;
    tord.reserve(n);    
    std::valarray<color_t> colors(boost::white_color, n);
    
    vertex_ref_map<color_t> cmap = &(colors[0]);
    
    // main
    
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
    
    // output
    
    if (is_dag)
    {
        size_t nt = tord.size();
        marray mRes = create_marray<int32_t>(1, nt);
        int32_t *res = mRes.data<int32_t>();
        
        for (size_t i = 1; i <= nt; ++i)
        {
            res[i-1] = tord[nt-i].index + 1;
        }
        
        plhs[0] = mRes.mx_ptr();
    }
    else
    {
        plhs[0] = create_mscalar<double>(-1).mx_ptr();
    }    
    
}


BCSMEX_MAINDEF
        
        

