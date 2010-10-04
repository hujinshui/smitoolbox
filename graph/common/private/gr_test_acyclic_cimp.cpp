/********************************************************************
 *
 *  gr_test_acyclic_cimp.cpp
 *
 *  The C++ mex implementation of gr_test_acyclic
 *
 *  Created by Dahua Lin, on Oct 3, 2010
 *
 *******************************************************************/

#include "../../clib/mgraph.h"
#include "../../clib/graph_search.h"

using namespace smi;


mxArray* revlist_to_matlab(const SeqList<int>& s)
{
    int n = s.size();
    mxArray *mx = create_matlab_matrix<int>(1, n);
    int *dst = (int*)mxGetData(mx);
    
    for (int i = 0; i < n; ++i)
    {
        dst[i] = s[n-(i+1)] + 1;
    }
    return mx;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mG(prhs[0]);
    
    RefGraph G = to_refgraph(mG);
    AdjList adjList(G);

    // main
    
    // do DFS
    
    DFSIteratorEx dfs_it(adjList);
    
    if (test_acyclic(G.nnodes(), dfs_it))
    {
        plhs[0] = create_matlab_scalar(true);
        if (nlhs >= 2)
        {
            plhs[1] = revlist_to_matlab(dfs_it.finish_order());
        }        
    }
    else
    {
        plhs[0] = create_matlab_scalar(false);
        if (nlhs >= 2)
            plhs[1] = create_empty_matrix();
    }
    
}


