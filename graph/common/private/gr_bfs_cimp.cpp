/********************************************************************
 *
 *  gr_bfs_cimp.cpp
 *
 *  The C++ mex implementation of gr_bfs_cimp
 *
 *  Created by Dahua Lin, on Oct 3, 2010
 *
 *******************************************************************/


#include "../../clib/mgraph.h"
#include "../../clib/graph_search.h"

using namespace smi; 
    
void do_bfs(const AdjList& adjList, mxArray*& mxVs)
{
    int *vs = new int[adjList.nnodes()];
    
    BFSIterator bfs_it(adjList);
    
    int v;
    int c = 0;
    while ((v = bfs_it.next()) >= 0)
    {
        vs[c++] = v;
    }
    
    mxVs = gindices_mrow(c, vs);
    delete[] vs;
}
    
    
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const MArray mG(prhs[0]);
    const MArray mSeeds(prhs[1]);
    
    // make graph
    
    RefGraph G = to_refgraph(mG);
    AdjList adjList(G);
    
    int n = mSeeds.nelems();
    const int *seeds = mSeeds.get_data<int>();
    
    // main
    
    if (nlhs <= 1)
    {
        do_bfs(adjList, plhs[0]);
    }    
}
    
