/********************************************************************
 *
 *  gr_traverse_cimp.cpp
 *
 *  The C++ mex implementation of gr_traverse.m
 *
 *  Created by Dahua Lin, on Oct 2, 2010
 *
 ********************************************************************/

#include "../../clib/mgraph_classic.h"

#include <string.h>

using namespace smi;

mxArray* do_bf_traverse(const GNeighborHood& G, int ns, const int *s)
{
    int n = G.nnodes();
    int *r = new int[n];
    bool *visited = new bool[n];
    ::memset(visited, 0, sizeof(bool) * n);
    
    int nr = 0;
    for (int i = 0; i < ns; ++i)
    {
        int v0 = s[i]-1;
        if (!visited[v0])
            nr += breadth_first_traverse(G, v0, r + nr, visited);
    }        
    mxArray *mxR = gindices_mrow(nr, r);
    
    delete[] r;
    delete[] visited;
    
    return mxR;
}

mxArray* do_df_traverse(const GNeighborHood& G, int ns, const int *s)
{
    int n = G.nnodes();
    int *r = new int[n];
    bool *visited = new bool[n];
    ::memset(visited, 0, sizeof(bool) * n);
    
    int nr = 0;
    for (int i = 0; i < ns; ++i)
    {
        int v0 = s[i]-1;
        if (!visited[v0])
            nr += depth_first_traverse(G, v0, r + nr, visited);
    }    
    mxArray *mxR = gindices_mrow(nr, r);
    
    delete[] r;
    delete[] visited;
    
    return mxR;
}

// main entry:
// Inputs:
//    [0] G:  the mgraph struct
//    [1] s:  starting node index array (one-based double)
//    [2] op: the opcode ('b' or 'd')
// Outpus:
//    [0] r:  the vector of traversed nodes (one-based)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxG = prhs[0];
    const mxArray *mxS = prhs[1];
    const mxArray *mxOp = prhs[2];
            
    MGraph G0(mxG);
    GNeighborHood G(G0, gnb_out());
        
    const int *s = (const int*)mxGetData(mxS);
    int ns = mxGetNumberOfElements(mxS);
    
    char op = char(*((const mxChar*)mxGetData(mxOp)));
    
    if (op == 'b')
    {
        plhs[0] = do_bf_traverse(G, ns, s);
    }
    else
    {
        plhs[0] = do_df_traverse(G, ns, s);
    }    
}




