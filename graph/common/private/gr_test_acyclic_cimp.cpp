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


#include <vector>
#include <functional>
#include <iterator>


using namespace smi;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mG(prhs[0]);
    
    RefGraph G = to_refgraph(mG);
    AdjList adjList(G);

    // main
    
    std::vector<int> vs;
    vs.reserve(G.nnodes());    
    
    DFSIterator dfs_it(adjList);
    
    if (test_acyclic(dfs_it, std::back_inserter(vs)))
    {
        plhs[0] = create_matlab_scalar(true);
        if (nlhs >= 2)
        {
            plhs[1] = iter_to_matlab_row(vs.rbegin(), (int)vs.size(), 
                    std::bind2nd(std::plus<int>(), 1));
        }        
    }
    else
    {
        plhs[0] = create_matlab_scalar(false);
        if (nlhs >= 2)
            plhs[1] = create_empty_matrix();
    }
    
}


