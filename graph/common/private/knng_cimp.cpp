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
#include <bcslib/graph/bgl_port.h>
#include <bcslib/matlab/mgraph.h>

#include <vector>
#include <functional>

using namespace bcs;
using namespace bcs::matlab;

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
gr_size_t knng_prune(const gr_wadjlist<TWeight, gr_undirected>& g, int K, 
        std::vector<int>& ss,
        std::vector<int>& st, 
        std::vector<TWeight>& sw)
{
    gr_size_t ns = 0;
    gr_size_t n = num_vertices(g);
    
    for (gr_size_t i = 0; i < n; ++i)
    {
        vertex_t v(i);
        
        gr_size_t d = out_degree(v, g);
        std::vector<ewinfo<TWeight> > tmp;
        tmp.reserve(d);        
        
        typedef typename gr_wadjlist<TWeight, gr_undirected>::adj_edge_iterator eiter_t;
        eiter_t eit;
        eiter_t eit_end;
        
        // collect local info
        
        for (rbind(eit, eit_end) = out_edges(v, g); eit != eit_end; ++eit) 
        {
            edge_t e = *eit;
            
            ewinfo<TWeight> ew;
            ew.s = source(e, g).index;
            ew.t = target(e, g).index;
            ew.w = g.weight_of(e);
            
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
inline void do_knng(const_mgraph mG, int K, bool use_min, int nlhs, mxArray *plhs[])
{
    gr_wadjlist<TWeight, gr_undirected> g = to_gr_wadjlist<TWeight, gr_undirected>(mG);
    
    gr_size_t n = num_vertices(g);
    
    std::vector<int32_t> ss;
    std::vector<int32_t> st;
    std::vector<TWeight> sw;
        
    ss.reserve(n * K);
    st.reserve(n * K);
    sw.reserve(n * K);
    
    gr_size_t ns = 0;
    
    if (use_min)
    {
        ns = knng_prune<TWeight, std::less>(g, K, ss, st, sw);
    }
    else
    {
        ns = knng_prune<TWeight, std::greater>(g, K, ss, st, sw);
    }
    
    marray mSS = create_marray<int32_t>(ns, 1);
    marray mST = create_marray<int32_t>(ns, 1);
    marray mSW = create_marray<TWeight>(ns, 1);
    
    int32_t *p_ss = mSS.data<int32_t>();
    int32_t *p_st = mST.data<int32_t>();
    TWeight *p_sw = mSW.data<TWeight>();
    
    for (gr_size_t i = 0; i < ns; ++i)
    {
        p_ss[i] = ss[i] + 1;
        p_st[i] = st[i] + 1;
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
 * [0]: g:  the gr_adjlist object
 * [1]: K:  the maximum number of neighbors 
 * [2]: use_min:  logical
 *
 * Out
 * [0]: ss: the sources of selected edges
 * [1]: st: the targets of selected edges
 * [2]: sw: the weights of selected edges
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const_mgraph mG(prhs[0]);    
    int K = (int)(const_marray(prhs[1]).get_scalar<double>());
    bool use_min = const_marray(prhs[2]).get_scalar<bool>();
    
    switch (mG.weight_type())
    {
        case mxDOUBLE_CLASS:
            do_knng<double>(mG, K, use_min, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            do_knng<float>(mG, K, use_min, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            do_knng<int32_t>(mG, K, use_min, nlhs, plhs);
            break;
            
        case mxUINT32_CLASS:
            do_knng<uint32_t>(mG, K, use_min, nlhs, plhs);
            break;
            
        default:
            throw mexception("knng_cimp:invalidarg", 
                "The weight value type should be double, single, int32, or uint32.");
    }
    
}

BCSMEX_MAINDEF






