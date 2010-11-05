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
#include "../../../graph/clib/graph_fvset.h"

#include <vector>

using namespace smi;

template<typename TWeight>
struct FVSVisitor
{
    FVSVisitor(graph_size_t n)
    {
        vs.reserve(n);
        scores.reserve(n);
    }    
    
    std::vector<vertex_t> vs;
    std::vector<TWeight> scores;
    
    void select_vertex(vertex_t v, TWeight s)
    {
        vs.push_back(v);
        scores.push_back(s);
    }
};



template<typename TWeight>
void do_gafvs(const matlab_graph_repr& gr, graph_size_t nmax, int nlhs, mxArray *plhs[])
{
    typedef CRefAdjList<TWeight, boost::undirected_tag> graph_t;
    graph_t g = gr.to_cref_wadjlist_ud<TWeight>();
    
    FVSVisitor<TWeight> vis(nmax);
    
    select_feedback_vertex(g, fvs_deg_computer<graph_t>(), nmax, vis);    
    
    plhs[0] = iter_to_matlab_row(vis.vs.begin(), vis.vs.size(), vertex_to_mindex());
    
    if (nlhs > 1)
    {
        plhs[1] = iter_to_matlab_row(vis.scores.begin(), vis.scores.size());
    }    
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
            
    switch (gr.weight_class())
    {
        case mxDOUBLE_CLASS:
            do_gafvs<double>(gr, nmax, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            do_gafvs<float>(gr, nmax, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            do_gafvs<int>(gr, nmax, nlhs, plhs);
            break;
            
        default:
            mexErrMsgIdAndTxt("gafva:invalidarg",
                    "The edge weights should be double, single, or int32.");
    }    
}


