/**********************************************************
 *
 *  mincut_bk_mex.cpp
 *
 *  The mex wrapper of Boykov's Max-flow Min-cut algorithm
 *
 *  Created by Dahua Lin, on Sep 26, 2010
 *
 **********************************************************/

#include <mex.h>

#include "maxflow-v3.01/graph.cpp"
#include "maxflow-v3.01/maxflow.cpp"


template<typename T>
Graph<T,T,T>* create_graph(int n, int m)
{    
    int nmax = n + 2;
    int mmax = 2 * (n + m + 1);
    
    Graph<T,T,T> *pG = new Graph<T,T,T>(nmax, mmax);
    pG->add_node(n);
    
    return pG;
}


template<typename T>
void release_graph(Graph<T,T,T> *pG)
{
    delete pG;
}


template<typename T>
struct Pair
{
    T first;
    T second;
};


template<typename T>
void add_tlinks(Graph<T,T,T> *pG, int n, const mxArray *mxTWs)
{    
    const Pair<T> *tws = (const Pair<T>*)mxGetData(mxTWs);
    
    for (int i = 0; i < n; ++i)
    {
        T w1 = tws[i].first;
        T w2 = tws[i].second;
        
        if (w1 != 0 && w2 != 0)
        {            
            pG->add_tweights(i, w1, w2);
        }                
    }
}

template<typename T>
void add_nlinks(Graph<T,T,T> *pG, int m, const mxArray *mxI, const mxArray *mxJ, const mxArray *mxV)
{
    const int *I = (const int*)mxGetData(mxI);
    const int *J = (const int*)mxGetData(mxJ);
    const T *V = (const T*)mxGetData(mxV);
    
    for (int k = 0; k < m; ++k)
    {
        int i = I[k];
        int j = J[k];
        T w = V[k];
        
        pG->add_edge(i-1, j-1, w, w);
    }        
}


template<typename T>
double do_mincut(int n, int m, 
        const mxArray *mxI, const mxArray *mxJ, const mxArray *mxV, 
        const mxArray *mxTWs, bool results[])
{
    // construct graph
    
    Graph<T,T,T> *pG = create_graph<T>(n, m);
    
    add_tlinks<T>(pG, n, mxTWs);
    add_nlinks<T>(pG, m, mxI, mxJ, mxV);
        
    // solve
    
    double mfv = (double)(pG->maxflow());
    
    // extract results
    
    for (int i = 0; i < n; ++i)
    {
        results[i] = (pG->what_segment(i) == Graph<T,T,T>::SINK);
    }        
    
    // release graph
    
    release_graph(pG);
    
    return mfv;
}


/************
 *
 *  main entry:
 *
 *  Input:
 *    [0,1,2]: I,J,V:   the nlinks
 *    [3] tws:      the weights to terminals (2 x n)
 *     
 *  Output
 *    [0] results:  the cut results (1 x n bool vector: 0 - source, 1 - sink)
 *    [1] maxflow:  the value of maximum flow (minimum cut)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{
    // take input
    
    const mxArray *mxI = prhs[0];
    const mxArray *mxJ = prhs[1];
    const mxArray *mxV = prhs[2];
    const mxArray *mxTWs = prhs[3];
    
    // size info
    
    int m = mxGetNumberOfElements(mxI);
    int n = mxGetN(mxTWs);
    mxClassID cid = mxGetClassID(mxV);
    
    double mfv = 0;
    mxArray *mxResults = mxCreateLogicalMatrix(1, n);
    bool *results = (bool*)mxGetData(mxResults);
        
    switch (cid)
    {
        case mxDOUBLE_CLASS:
            mfv = do_mincut<double>(n, m, mxI, mxJ, mxV, mxTWs, results);
            break;
        case mxSINGLE_CLASS:
            mfv = do_mincut<float>(n, m, mxI, mxJ, mxV, mxTWs, results);
            break;
        case mxINT32_CLASS:
            mfv = do_mincut<int>(n, m, mxI, mxJ, mxV, mxTWs, results);
            break;            
    }    
        
    plhs[0] = mxResults;
    if (nlhs >= 2)
    {
        plhs[1] = mxCreateDoubleScalar(mfv);
    }
}









