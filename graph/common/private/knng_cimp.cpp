/**********************************************************
 *
 *  knng_cimp.cpp
 *
 *  The C++ mex implementation for part of knng.m
 *
 *  Created by Dahua Lin, on Nov 13, 2010
 *
 **********************************************************/ 

#include <bcslib/matlab/bcs_mex.h>
#include "../../clib/smi_graph_mex.h"

#include <vector>
#include <functional>

using namespace bcs;
using namespace bcs::matlab;
using namespace smi;

template<typename T>
struct ewinfo
{
    int32_t s;
    int32_t t;
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
gindex_t knng_prune(const DGraph& g, 
        const TWeight *weights,
        int K, 
        std::vector<int>& ss,
        std::vector<int>& st, 
        std::vector<TWeight>& sw)
{
    gindex_t ns = 0;
    gindex_t n = g.nnodes();
    
    for (gindex_t v = 1; v <= n; ++v)
    {        
        gindex_t d = g.out_degree(v);
        std::vector<ewinfo<TWeight> > tmp;
        tmp.reserve(d);        
                
        // collect local info
        
        const gindex_t *E = g.out_edges(v);
        
        for (gindex_t j = 0; j < d; ++j) 
        {
            gindex_t eidx = E[j];
            const DEdge& e = g.edge(eidx);
            
            ewinfo<TWeight> ew;
            ew.s = e.s;
            ew.t = e.t;
            ew.w = weights[eidx];
            
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
        
        ns += k;        
    }
    
    return ns;
}


template<typename TWeight>
inline void do_knng(const DGraph& g, const TWeight *weights, 
        int K, bool use_min, int nlhs, mxArray *plhs[])
{    
    gindex_t n = g.nnodes();
    
    std::vector<int32_t> ss;
    std::vector<int32_t> st;
    std::vector<TWeight> sw;
        
    ss.reserve(n * K);
    st.reserve(n * K);
    sw.reserve(n * K);
    
    gindex_t ns = 0;
    
    if (use_min)
    {
        ns = knng_prune<TWeight, std::less>(g, weights, K, ss, st, sw);
    }
    else
    {
        ns = knng_prune<TWeight, std::greater>(g, weights, K, ss, st, sw);
    }
    
    marray mSS = create_marray<int32_t>(ns, 1);
    marray mST = create_marray<int32_t>(ns, 1);
    marray mSW = create_marray<TWeight>(ns, 1);
    
    int32_t *p_ss = mSS.data<int32_t>();
    int32_t *p_st = mST.data<int32_t>();
    TWeight *p_sw = mSW.data<TWeight>();
    
    for (gindex_t i = 0; i < ns; ++i)
    {
        p_ss[i] = ss[i];
        p_st[i] = st[i];
        p_sw[i] = sw[i];
    }
        
    plhs[0] = mSS.mx_ptr();    
    plhs[1] = mST.mx_ptr();    
    plhs[2] = mSW.mx_ptr();
}



/**
 * main entry:
 *
 * In
 * [0]: g:  the graph struct
 * [1]: w:  the edge weights
 * [2]: K:  the maximum number of neighbors 
 * [3]: use_min:  logical
 *
 * Out
 * [0]: ss: the sources of selected edges
 * [1]: st: the targets of selected edges
 * [2]: sw: the weights of selected edges
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const_marray mG(prhs[0]);
    DGraph g(get_dgraph_spec(mG));    
    
    const_marray mW(prhs[1]);
    
    int K = (int)(const_marray(prhs[2]).get_scalar<double>());
    bool use_min = const_marray(prhs[3]).get_scalar<bool>();
    
    switch (mW.class_id())
    {
        case mxDOUBLE_CLASS:
            do_knng<double>(g, mW.data<double>(), K, use_min, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            do_knng<float>(g, mW.data<float>(), K, use_min, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            do_knng<int32_t>(g, mW.data<int32_t>(), K, use_min, nlhs, plhs);
            break;
            
        case mxUINT32_CLASS:
            do_knng<uint32_t>(g, mW.data<uint32_t>(), K, use_min, nlhs, plhs);
            break;
            
        default:
            throw mexception("knng_cimp:invalidarg", 
                "The weight value type should be double, single, int32, or uint32.");
    }
    
}

BCSMEX_MAINDEF






