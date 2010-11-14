/**********************************************************
 *
 *  knng_cimp.cpp
 *
 *  The C++ mex implementation for part of knng.m
 *
 *  Created by Dahua Lin, on Nov 13, 2010
 *
 **********************************************************/ 

#include <mex.h>

#include "../../clib/graph_mex.h"
#include <vector>
#include <functional>

using namespace smi;

template<typename T>
struct ewinfo
{
    int s;
    int t;
    T w;
    
    bool operator < (const ewinfo<T>& rhs) const
    {
        return w < rhs.w;
    }
    
    bool operator == (const ewinfo<T>& rhs) const
    {
        return w == rhs.w;
    }
    
    bool operator > (const ewinfo<T>& rhs) const
    {
        return w > rhs.w;
    }
};



template<typename TWeight, template<typename U> class TComp>
void knng_prune(const CRefAdjList<TWeight>& g, int K, 
        std::vector<int>& ss,
        std::vector<int>& st, 
        std::vector<TWeight>& sw)
{
    graph_size_t n = num_vertices(g);
    
    for (graph_size_t i = 0; i < n; ++i)
    {
        vertex_t v(i);
        
        graph_size_t d = out_degree(v, g);
        std::vector<ewinfo<TWeight> > tmp;
        tmp.reserve(d);        
        
        typedef typename CRefAdjList<TWeight>::out_edge_iterator eiter_t;
        eiter_t eit;
        eiter_t eit_end;
        
        // collect local info
        
        for (boost::tie(eit, eit_end) = out_edges(v, g); eit != eit_end; ++eit) 
        {
            edge_t e = *eit;
            
            ewinfo<TWeight> ew;
            ew.s = source(e, g).i;
            ew.t = target(e, g).i;
            ew.w = g.get_weight(e);
            
            tmp.push_back(ew);
        }
        
        // select
        
        int k = K;
        
        if (K < (int)tmp.size()) 
        {
            std::nth_element(tmp.begin(), tmp.begin() + k, tmp.end(), 
                    TComp<ewinfo<TWeight> >() );
        }
        else {
            k = (int)tmp.size();
        }
        
        // store
        
        for (int i = 0; i < k; ++i) 
        {
            const ewinfo<TWeight>& ew = tmp[i];
            
            ss.push_back(ew.s);
            st.push_back(ew.t);
            sw.push_back(ew.w);            
        } 
    }
}


template<typename TWeight>
inline void main_delegate(const matlab_graph_repr& gr, int K, bool use_min, 
        int nlhs, mxArray *plhs[])
{
    CRefAdjList<TWeight> g = gr.to_cref_wadjlist<TWeight>();
    
    graph_size_t n = num_vertices(g);
    
    std::vector<int> ss;
    std::vector<int> st;
    std::vector<TWeight> sw;
        
    ss.reserve(n * K);
    st.reserve(n * K);
    sw.reserve(n * K);
    
    if (use_min)
    {
        knng_prune<TWeight, std::less>(g, K, ss, st, sw);
    }
    else
    {
        knng_prune<TWeight, std::greater>(g, K, ss, st, sw);
    }
        
    plhs[0] = iter_to_matlab_column(ss.begin(), ss.size(), 
            std::bind2nd(std::plus<int>(), 1));
    
    plhs[1] = iter_to_matlab_column(st.begin(), st.size(), 
            std::bind2nd(std::plus<int>(), 1));
    
    plhs[2] = iter_to_matlab_column(sw.begin(), sw.size());
}



/**
 * main entry:
 *
 * In
 * [0]: g:  the gr_adjlist object
 * [1]: K:  the maximum number of neighbors 
 * [2]: use_min:  logical
 *
 * Out
 * [0]: ss: the sources of selected edges
 * [1]: st: the targets of selected edges
 * [2]: sw: the weights of selected edges
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    matlab_graph_repr gr(prhs[0]);    
    MArray mK(prhs[1]);
    MArray mUMin(prhs[2]);
    
    int K = (int)mK.get_scalar<double>();
    bool use_min = mUMin.get_scalar<bool>();
    
    switch (gr.weight_class())
    {
        case mxDOUBLE_CLASS:
            main_delegate<double>(gr, K, use_min, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            main_delegate<float>(gr, K, use_min, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            main_delegate<int>(gr, K, use_min, nlhs, plhs);
            break;
            
        case mxUINT32_CLASS:
            main_delegate<unsigned int>(gr, K, use_min, nlhs, plhs);
            break;
            
        default:
            mexErrMsgIdAndTxt("knng_cimp:invalidarg", 
                    "The weight value should be double, single, int32, or uint32.");
    }
    
}








