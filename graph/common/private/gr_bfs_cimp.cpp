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
    
void do_bfs(const AdjList& adjList, int ns, const int *seeds, mxArray*& mxVs)
{
    BFSIterator bfs_it(adjList);
    
    for (int i = 0; i < ns; ++i) 
        bfs_it.add_seed(seeds[i]);
    
    SeqList<int> vs(adjList.nnodes());
    
    int v = -1;
    while ((v = bfs_it.next()) >= 0)
    {
        vs.add(v);
    }
    
    mxVs = gindices_mrow(vs.size(), vs.data());    
}
        
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const MArray mG(prhs[0]);
    const MArray mSeeds(prhs[1]);
    
    // make graph
    
    RefGraph G = to_refgraph(mG);
    AdjList adjList(G);
    
    int ns = mSeeds.nelems();
    const int *seeds = mSeeds.get_data<int>();
    
    // main
    
    if (nlhs <= 1)
    {
        do_bfs(adjList, ns, seeds, plhs[0]);
    }    
}
    
