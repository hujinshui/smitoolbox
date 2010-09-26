/********************************************************************
 *
 *  Wrapper of Boycov's implementation of Maximum-flow Graph cut
 *
 *  Created by Dahua Lin, on Oct 11, 2009
 *  Modified by Dahua Lin, on Nov 6, 2009 (change weight type to double)
 *
 ********************************************************************/


#include <mex.h>

// These .cpp files are actually definitions of template classes and 
// template functions, so these files should be included in code but not
// linked as libraries

// #include <graph.h>
#include "maxflow-v3/graph.cpp"
#include "maxflow-v3/maxflow.cpp"


typedef double weight_type;

typedef Graph<weight_type, weight_type, weight_type> WGraph;


// The struct to represent an edge
struct Edge
{
    int i1;
    int i2;
};


/**
 * The function to create the graph from edges
 *
 * Remarks:
 *  - the terminal weight: positive --> close to source, negative --> close to sink 
 *  - all edges are symmetry, cap = rev_cap = neighbor weight
 */
WGraph *create_wgraph(
        int n,              /* the number of nodes */
        int m,              /* the number of edges */
        const Edge *edges,  /* the array of edges */
        const weight_type* tweights,    /* the array of terminal weights (n) */ 
        const weight_type* nweights     /* the array of neighbor weights (m) */
        ) 
{
    int nmax = n + 2;
    int mmax = m + n * 2;
    
    WGraph *pG = new WGraph(nmax, mmax);
    
    // add notes and terminal links
    
    pG->add_node(n);
    for (int i = 0; i < n; ++i)
    {
        weight_type tw = tweights[i];
        if (tw > 0) 
            pG->add_tweights(i, tw, 0);
        else if (tw < 0)
            pG->add_tweights(i, 0, -tw);
    }
    
    // add edges (neighboring links)
    for (int i = 0; i < m; ++i)
    {
        const Edge& e = edges[i];
        weight_type nw = nweights[i];
        
        pG->add_edge(e.i1 - 1, e.i2 - 1, nw, nw);
    }    
    
    return pG;
}


/**
 * Extract results from the solution
 * 
 * In the results: true -> source, false -> sink
 */
void extract_results(int n, WGraph *pG, bool *results)
{
    for (int i = 0; i < n; ++i)
    {
        results[i] = (pG->what_segment(i) == WGraph::SINK);
    }
}


/**
 * The main entry
 *
 * Input 
 *  - n:        the number of nodes [double scalar]
 *  - edges:    the array of edges (using one-based indices of nodes) [2 x m int32 matrix]
 *  - tweights: the array of terminal weights [1 x n double vector]
 *  - nweights: the array of neighbor weights [1 x m double vector]
 *
 * Output
 *  - results:  the result vector [1 x n logical array: 1 -> source, 0 -> sink]
 *  - maxflow:  the maximum flow (minimum cut)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
        
    const mxArray *mxN = prhs[0];
    const mxArray *mxEdges = prhs[1];
    const mxArray *mxTweights = prhs[2];
    const mxArray *mxNweights = prhs[3];
    
    // verify input arguments
        
    int n = (int)mxGetScalar(mxN);
    int m = mxGetN(mxEdges);
    
    const Edge *edges = (const Edge*)mxGetData(mxEdges);
    const weight_type *tweights = (const weight_type*)mxGetData(mxTweights);
    const weight_type *nweights = (const weight_type*)mxGetData(mxNweights);
    
    // prepare output
    
    mxArray *mxResults = mxCreateLogicalMatrix(1, n);
    bool *results = (bool*)mxGetData(mxResults);
    
    // create graph and solve
    
    WGraph *pG = create_wgraph(n, m, edges, tweights, nweights);
    weight_type mfv = pG->maxflow();
    
    extract_results(n, pG, results);
    
    // output
    
    delete pG;
    plhs[0] = mxResults;    
    
    if (nlhs >= 2)
    {
        plhs[1] = mxCreateDoubleScalar(mfv);
    }
}











