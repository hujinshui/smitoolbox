/**********************************************************
 *
 *  kolmogorov_mincut_cimp.cpp
 *
 *  A C++ mex wrapper of Kolmogorov's mincut implementation
 *  (based on maxflow v3.01)
 *
 * Created by Dahua Lin, on Jan 29, 2012
 *
 **********************************************************/

#include <mex.h>
#include "kolmogorov_maxflow.h"

inline mxClassID decide_capacity_type(
        const mxArray *mxNc, 
        const mxArray *mxRc, 
        const mxArray *mxTcSrc,
        const mxArray *mxTcSnk)
{
    mxClassID cid_nc = mxGetClassID(mxNc);
    mxClassID cid_rc = mxGetClassID(mxRc);
    mxClassID cid_tsrc = mxGetClassID(mxTcSrc);
    mxClassID cid_tsnk = mxGetClassID(mxTcSnk);
    
    if (!(cid_nc == cid_rc && cid_rc == cid_tsrc && cid_tsrc == cid_tsnk))
    {
        mexErrMsgIdAndTxt("kolmogorov_mincut:invalidarg", 
            "The neighbor-link capacities and terminal-link capacities have different types.");
    }
    
    return cid_nc;
}


template<typename T>
void do_cut(mxClassID cid, int n, int m, const int *sv, const int *tv, 
        const T* nb_cap, const T* rv_cap, const T* tsrc_cap, const T* tsnk_cap, 
        mxArray *plhs[])
{
    typedef vkol::Graph<T,T,T> vk_graph_t;
    
    // construct graph
    
    vk_graph_t G(n, m);    
    G.add_node(n);
    
    // add terminal capacities
    
    for (int i = 0; i < n; ++i)
    {
        T src_c = tsrc_cap[i];
        T snk_c = tsnk_cap[i];
        
        if (src_c > 0 || snk_c > 0)
        {
            G.add_tweights(i, src_c, snk_c);
        }
    }
    
    // add neighboring capacities
    
    for (int i = 0; i < m; ++i)
    {
        T nb_c = nb_cap[i];
        T rv_c = rv_cap[i];
        
        if (nb_c > 0 || rv_c > 0)
        {
            G.add_edge(sv[i]-1, tv[i]-1, nb_c, rv_c); 
        }
    }
    
    // solve max-flow
    
    T fv = G.maxflow();
    
    // extract results
    
    mxArray *mxL = mxCreateNumericMatrix(1, n, mxINT32_CLASS, mxREAL);
    int *L = (int*)mxGetData(mxL);
    
    for (int i = 0; i < n; ++i)
    {
        L[i] = G.what_segment_ex(i);
    }
        
    mxArray *mxFlow = mxCreateNumericMatrix(1, 1, cid, mxREAL);    
    *((T*)mxGetData(mxFlow)) = fv;
    
    plhs[0] = mxL;
    plhs[1] = mxFlow;
}



/**
 * Input 
 *  [0]  n:     The number of vertices [double]
 *  [1]  m:     The number of links [double]
 *  [2]  sv:    The vector of source vertices of edges [int32 one-based]
 *  [3]  tv:    The vector of target vertices of edges [int32 one-based]
 *  [4]  nc:    The capacities of neighboring-links
 *  [5]  rc:    The reverse capacities of neighboring-links
 *  [6]  tc_src:  The capacities of links from source
 *  [7]  tc_snk:  The capacities of links from sink
 *  
 *  Note the type of nc can be double|single|int32|uint32, tc should be
 *  of the same type as nc 
 *
 * Output
 *  [0]  L:     The labeling of the nodes [int32 1 x n]
 *                  - 1:    assigned to source
 *                  - -1:   assigned to sink
 *                  - 0:    can be assigned to either source or sink
 *
 *  [1]  maxflow:   The value of maxflow           
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxN = prhs[0];
    const mxArray *mxM = prhs[1];
    const mxArray *mxSv = prhs[2];
    const mxArray *mxTv = prhs[3];
    
    const mxArray *mxNc = prhs[4];
    const mxArray *mxRc = prhs[5];
    const mxArray *mxTcSrc = prhs[6];
    const mxArray *mxTcSnk = prhs[7];
    
    int n = (int)mxGetScalar(mxN);
    int m = (int)mxGetScalar(mxM);
    const int *sv = (const int*)mxGetData(mxSv);
    const int *tv = (const int*)mxGetData(mxTv);
        
    mxClassID cid = decide_capacity_type(mxNc, mxRc, mxTcSrc, mxTcSnk);
    
    switch (cid)
    {
        case mxDOUBLE_CLASS:
            {
                const double *nb_cap = (const double*)mxGetData(mxNc);
                const double *rv_cap = (const double*)mxGetData(mxRc);
                const double *tsrc_cap = (const double*)mxGetData(mxTcSrc);
                const double *tsnk_cap = (const double*)mxGetData(mxTcSnk);            
                do_cut(cid, n, m, sv, tv, nb_cap, rv_cap, tsrc_cap, tsnk_cap, plhs);            
            }
            break;
            
        case mxSINGLE_CLASS:                        
            {
                const float *nb_cap = (const float*)mxGetData(mxNc);
                const float *rv_cap = (const float*)mxGetData(mxRc);
                const float *tsrc_cap = (const float*)mxGetData(mxTcSrc);
                const float *tsnk_cap = (const float*)mxGetData(mxTcSnk);   
                do_cut(cid, n, m, sv, tv, nb_cap, rv_cap, tsrc_cap, tsnk_cap, plhs);  
            }
            break;
            
        case mxINT32_CLASS:                        
            {
                const int *nb_cap = (const int*)mxGetData(mxNc);
                const int *rv_cap = (const int*)mxGetData(mxRc);
                const int *tsrc_cap = (const int*)mxGetData(mxTcSrc);
                const int *tsnk_cap = (const int*)mxGetData(mxTcSnk);   
                do_cut(cid, n, m, sv, tv, nb_cap, rv_cap, tsrc_cap, tsnk_cap, plhs);  
            }
            break;
            
        case mxUINT32_CLASS:  
            {
                const unsigned int *nb_cap = (const unsigned int*)mxGetData(mxNc);
                const unsigned int *rv_cap = (const unsigned int*)mxGetData(mxRc);
                const unsigned int *tsrc_cap = (const unsigned int*)mxGetData(mxTcSrc);
                const unsigned int *tsnk_cap = (const unsigned int*)mxGetData(mxTcSnk);   
                do_cut(cid, n, m, sv, tv, nb_cap, rv_cap, tsrc_cap, tsnk_cap, plhs); 
            }
            break;
            
        default:
            mexErrMsgIdAndTxt("kolmogorov_mincut:invalidarg", 
                "The capacity type is unsupported.");
    }
    
}

