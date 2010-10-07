/********************************************************************
 *
 *  gr_conncomp_cimp.cpp
 *
 *  The C++ mex implementation of gr_conncomp
 *
 *  Created by Dahua Lin, on Oct 4, 2010
 *
 ********************************************************************/


#include "../../clib/mgraph.h"
#include "../../clib/graph_search.h"

#include <vector>
#include <iterator>

using namespace smi;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mG(prhs[0]);
    RefGraph G = to_refgraph(mG);
    AdjList adjlist(G);
    
    // find components
    
    int n = G.nnodes();
    std::vector<int> csizs;
    std::vector<int> vs;
    vs.reserve(n);
    
    int c = get_connected_components(adjlist, 
            std::back_inserter(csizs), std::back_inserter(vs));
 
    // make output
    
    mxArray *mxCCs = mxCreateCellMatrix(1, c);
    const int *p = &(vs[0]);
    for (int i = 0; i < c; ++i)
    {
        mxArray *mxC = gindices_mrow(csizs[i], p);
        p += csizs[i];
        
        mxSetCell(mxCCs, i, mxC);
    }
    
    plhs[0] = mxCCs;    
}


