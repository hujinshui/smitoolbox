/********************************************************************
 *
 *  gafvs_cimp.cpp
 *
 *  The C++ mex implementation of gafvs.m
 *
 *  Created by Dahua Lin, on Nov 4, 2010
 *
 ********************************************************************/

#include "../../../graph/clib/graph_mex.h"

#include "feedback_vertex_set.h"

#include <vector>
#include <iterator>

using namespace smi;


template<typename TWeight>
void do_gafvs(const matlab_graph_repr& gr, graph_size_t nmax, std::vector<vertex_t>& fvs)
{
    typedef CRefAdjList<TWeight, boost::undirected_tag> graph_t;
    graph_t g = gr.to_cref_wadjlist_ud<TWeight>();
    
    select_feedback_vertex(g, fvs_deg_computer<graph_t>(), nmax, 
            std::back_inserter(fvs));    
}


/**
 * main entry:
 *
 * Input
 *   [0]: g: the graph (in form of gr_adjlist undirected and weighted)
 *   [1]: nmax: the maximum number of feedback vertices to extract (double)
 *
 * Output:
 *   [0]: fvs:  the row vector of extracted feedback vertices
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mG(prhs[0]);
    MArray mNmax(prhs[1]);
    
    matlab_graph_repr gr(mG);
    graph_size_t nmax = (graph_size_t)mNmax.get_double_scalar();
    
    std::vector<vertex_t> fvs;
        
    switch (gr.weight_class())
    {
        case mxDOUBLE_CLASS:
            do_gafvs<double>(gr, nmax, fvs);
            break;
            
        case mxSINGLE_CLASS:
            do_gafvs<float>(gr, nmax, fvs);
            break;
            
        case mxINT32_CLASS:
            do_gafvs<int>(gr, nmax, fvs);
            break;
            
        default:
            mexErrMsgIdAndTxt("gafva:invalidarg",
                    "The edge weights should be double, single, or int32.");
    }
    
    plhs[0] = iter_to_matlab_row(fvs.begin(), fvs.size(), vertex_to_mindex());
}


