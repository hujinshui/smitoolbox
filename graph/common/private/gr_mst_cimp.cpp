/********************************************************************
 *
 *  gr_mst_cimp.cpp
 *
 *  The C++ mex implementation for MST algorithms
 *
 *  Created by Dahua Lin, on Oct 8, 2010
 *
 ********************************************************************/

#include "../../clib/mgraph.h"
#include "../../clib/graph_mst.h"

#include <iterator>

using namespace smi;


template<typename T>
void create_edge_array(const std::vector<WEdge<T> >& edges, 
        mxArray*& mxI, mxArray*& mxJ, mxArray*& mxW)
{
    int n = (int)edges.size();
    
    mxI = create_matlab_matrix<int>(n,1);
    mxJ = create_matlab_matrix<int>(n,1);
    mxW = create_matlab_matrix<T>(n,1);
    
    int *I = (int*)mxGetData(mxI);
    int *J = (int*)mxGetData(mxJ);
    T *W = (T*)mxGetData(mxW);
    
    for (int i = 0; i < n; ++i)
    {
        const WEdge<T>& e = edges[i];
        
        I[i] = e.s + 1;
        J[i] = e.t + 1;
        W[i] = e.w;
    }    
}



template<typename T>
void do_prim_mst(const MArray& mG, int root, 
        mxArray*& mxI, mxArray*& mxJ, mxArray*& mxW)
{
    RefWGraph<T> G = to_refwgraph<T>(mG);
    WAdjList<T> adjList(G);
    
    std::vector<WEdge<T> > edges;
    edges.reserve(G.nnodes());
        
    Prim_MSTIterator<T> mst_it(adjList);
    
    if (root >= 0)
    {
        mst_it.initialize(root);
        mst_it.solve_all(std::back_inserter(edges));
    }
    else
    {
        int n = G.nnodes();
        for (int i = 0; i < n; ++i)
        {
            if (!mst_it.is_included(i))
            {
                mst_it.initialize(i);
                mst_it.solve_all(std::back_inserter(edges));
            }
        }
    }    
    
    create_edge_array(edges, mxI, mxJ, mxW);
}


template<typename T>
void do_kruskal_mst(const MArray& mG, 
        mxArray*& mxI, mxArray*& mxJ, mxArray*& mxW)
{
    RefWGraph<T> G = to_refwgraph<T>(mG);
    WAdjList<T> adjList(G);
    
    std::vector<WEdge<T> > edges;
    edges.reserve(G.nnodes());
        
    Kruskal_MSTIterator<T> mst_it(adjList);
    
    mst_it.initialize();
    mst_it.solve_all(std::back_inserter(edges));
    
    create_edge_array(edges, mxI, mxJ, mxW);
}



template<typename T>
void do_mst(const MArray& mG, char code, const mxArray *prhs[], 
        mxArray*& mxI, mxArray*& mxJ, mxArray*& mxW)
{
    if (code == 'p')
    {        
        MArray mRoot(prhs[2]);
        
        int root = (int)(mRoot.get_double_scalar()) - 1;        
        do_prim_mst<T>(mG, root, mxI, mxJ, mxW);
    }
    else if (code == 'k')
    {
        do_kruskal_mst<T>(mG, mxI, mxJ, mxW);
    }
}



// Main entry:
//  Input:
//      [0]: G:     the graph
//      [1]: code:  'p' -> Prim, 'k' -> Kruskal
//  Output:
//      [0]: I:     parent nodes
//      [1]: J:     target nodes
//      [2]: W:     weights
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    MArray mG(prhs[0]);    
    mxClassID wcid = get_graph_weight_class(mG);
    
    MArray mCode(prhs[1]);
    char code = (char)(mCode.get_scalar<mxChar>());
    
    switch (wcid)
    {
        case mxDOUBLE_CLASS:
            do_mst<double>(mG, code, prhs, plhs[0], plhs[1], plhs[2]);
            break;
            
        case mxSINGLE_CLASS:
            do_mst<float>(mG, code, prhs, plhs[0], plhs[1], plhs[2]);
            break;
            
        case mxINT32_CLASS:
            do_mst<int>(mG, code, prhs, plhs[0], plhs[1], plhs[2]);
            break;
            
        default:
            mexErrMsgIdAndTxt("gr_mst:invalidarg", 
                    "The weights should be either double, single, or int32.");
    }
}








