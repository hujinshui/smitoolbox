/**********************************************************
 *
 *  mincut_kolmogorov_cimp.cpp
 *
 *  The mex wrapper of Boykov's Max-flow Min-cut algorithm
 *
 *  Created by Dahua Lin, on Sep 26, 2010
 *
 **********************************************************/

#include "../../clib/graph_mex.h"

#include "../../clib/graph_mincut_kolmogorov.h"

using namespace smi;


template<typename TWeight>
void main_delegate(const matlab_graph_repr& gr, 
        const MArray& mSrcWs, const MArray& mSinkWs, 
        int nlhs, mxArray *plhs[])
{
    CRefEdgeList<TWeight, boost::undirected_tag> g = 
            gr.to_cref_wedgelist_ud<TWeight>();
    
    graph_size_t n = num_vertices(g);
    
    const TWeight *src_ws = mSrcWs.get_data<TWeight>();
    const TWeight *sink_ws = mSinkWs.get_data<TWeight>();
    
    mxArray *mxR = mxCreateLogicalMatrix(1, n);
    bool *results = (bool*)mxGetData(mxR);
    
    TWeight mfv = mincut_kolmogorov(g, src_ws, sink_ws, results);
    
    plhs[0] = mxR;
    if (nlhs >= 2)
        plhs[1] = create_matlab_scalar<TWeight>(mfv);
}




/************
 *
 *  main entry:
 *
 *  Input:
 *    [0] g:        the graph (in form of edgelist or adjlist, undirected)
 *    [1] src_ws:   the weights to source
 *    [2] sink_ws:  the weights to sink
 *     
 *  Output
 *    [0] results:  the cut results (1 x n bool vector: 0 - source, 1 - sink)
 *    [1] maxflow:  the value of maximum flow (minimum cut)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{
    // take input
    
    matlab_graph_repr gr(prhs[0]);    
    
    MArray mSrcWs(prhs[1]);
    MArray mSinkWs(prhs[2]);           
                
    // delegate
    
    switch (gr.weight_class())
    {
        case mxDOUBLE_CLASS:
            main_delegate<double>(gr, mSrcWs, mSinkWs, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            main_delegate<float>(gr, mSrcWs, mSinkWs, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            main_delegate<int>(gr, mSrcWs, mSinkWs, nlhs, plhs);
            break;
            
        case mxUINT32_CLASS:
            main_delegate<unsigned int>(gr, mSrcWs, mSinkWs, nlhs, plhs);
            break;
            
        default:
            mexErrMsgIdAndTxt("mincut_kolmogorov:invalidarg", 
                    "The weight value should be double, single, int32, or uint32.");
    }
}

