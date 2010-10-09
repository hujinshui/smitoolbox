/********************************************************************
 *
 *  gr_adjlist_cimp.cpp
 *
 *  The C++ mex implementation of adjacency list construction
 *
 *  Created by Dahua Lin, on Oct 8, 2010
 *
 *******************************************************************/


#include "../../clib/graph_mex.h"

#include <string.h>


using namespace smi;

template<typename TWeight>
bool has_weight() { return true; }

template<>
bool has_weight<no_edge_weight>() { return false; }        




// The workspace for adjlist construction
template<typename TWeight>
struct adjlist_ws
{
    graph_size_t n;
    graph_size_t m;
    
    vertex_t *s;
    vertex_t *t;
    TWeight  *w;
    
    graph_size_t *o_degs;
    graph_size_t *o_offsets;
    
    adjlist_ws() : n(0), m(0), s(0), t(0), w(0), o_degs(0), o_offsets(0) 
    {        
    }
       
    ~adjlist_ws()
    {
        if (s) delete[] s;
        if (t) delete[] t;
        if (w) delete[] w;
        
        if (o_degs) delete[] o_degs;
        if (o_offsets) delete[] o_offsets;
    } 
    
    
    template<class TGraph>
    void construct_from_edges(const TGraph& g)
    {
        typedef typename TGraph::edge_iterator edge_iterator;
        
        graph_size_t n = num_vertices(g);
        
        // scan edge numbers
        
        o_degs = new graph_size_t[n];
        o_offsets = new graph_size_t[n];
        
        ::memset(o_degs, 0, sizeof(graph_size_t) * n);
        
        edge_iterator it_begin, it_end;
        boost::tie(it_begin, it_end) = edges(g);
        
        for (edge_iterator it = it_begin; it != it_end; ++it)
        {
            ++ o_degs[source(*it, g).i];
        }
        
        graph_size_t m = 0;
        for (graph_size_t i = 0; i < n; ++i)
        {
            o_offsets[i] = m;
            m += o_degs[i];
        }
        
        // fill adjacency information
        
        s = new vertex_t[m];
        t = new vertex_t[m];                

        if (has_weight<TWeight>())
            w = new TWeight[m];        
        
        graph_size_t *c = new graph_size_t[n];
        ::memcpy(c, o_offsets, n * sizeof(graph_size_t));
        
        for (edge_iterator it = it_begin; it != it_end; ++it)
        {
            edge_t e = *it;
            vertex_t e_s = source(e, g);
            vertex_t e_t = target(e, g);
            
            graph_size_t idx = (c[e_s.i] ++);
                        
            s[idx] = e_s;
            t[idx] = e_t;
            if (w) w[idx] = g.get_weight(e);
        }
        
        delete[] c;
    }
       
    
    mxArray* output_matlab()
    {
        return create_matlab_graph_struct(n, m, s, t, w, o_degs, o_offsets);
    }
           
    
};



template<typename TWeight> 
mxArray* make_adjlist_e(const matlab_graph_repr& mg)
{    
    if (!has_weight<TWeight>())
    {
        CRefEdgeList<no_edge_weight> edges = mg.to_cref_edgelist(); 
        
        adjlist_ws<no_edge_weight> aws;
        aws.construct_from_edges(edges);
        
        return aws.output_matlab();
    }
    else
    {
        CRefEdgeList<TWeight> edges = mg.to_cref_wedgelist<TWeight>();
        
        adjlist_ws<TWeight> aws;
        aws.construct_from_edges(edges);
        
        return aws.output_matlab();
    }
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    matlab_graph_repr mg(prhs[0]);
    
    if (mg.has_weight())
    {
        switch (mg.weight_class())
        {
            case mxDOUBLE_CLASS:
                plhs[0] = make_adjlist_e<double>(mg);
                break;
                
            case mxSINGLE_CLASS:
                plhs[0] = make_adjlist_e<float>(mg);
                break;
                
            case mxINT32_CLASS:
                plhs[0] = make_adjlist_e<int>(mg);
                break;
                
            case mxUINT32_CLASS:
                plhs[0] = make_adjlist_e<unsigned int>(mg);
                break;
                
            default:
                mexErrMsgIdAndTxt("gr_adjlist:invalidarg", 
                        "The weights should be of class double, single, int32, or uint32.");
        }
    }
    else
    {
        plhs[0] = make_adjlist_e<no_edge_weight>(mg);
    }
}





