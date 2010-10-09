/********************************************************************
 *
 *  graph_mex.h
 *
 *  The mex interface for graph classes
 *
 *  Created by Dahua Lin, on Oct 8, 2010
 *
 ********************************************************************/


#ifndef SMI_GRAPH_MEX_H
#define SMI_GRAPH_MEX_H


#include "../../base/clib/matlab_types.h"
#include "../../base/clib/marray.h"

#include "graph_adjlist.h"


namespace smi
{

struct matlab_graph_repr
{
    const MArray& mG;
    
    matlab_graph_repr(const MArray& G) : mG(G) { }
    
    // field extraction
    
    bool has_weight() const
    {
        return !mG.get_field("w").is_empty();
    }
    
    mxClassID weight_class() const
    {
        return mG.get_field("w").class_id();
    }
    
    graph_size_t n() const
    {
        return mG.get_field("n").get_double_scalar();
    }
    
    graph_size_t m() const
    {
        return mG.get_field("m").get_double_scalar();
    }
    
    const vertex_t* s() const
    {
        return mG.get_field("s").get_data<vertex_t>();
    }
    
    const vertex_t *t() const
    {
        return mG.get_field("t").get_data<vertex_t>();
    }
    
    template<typename TWeight>
    const TWeight *w() const
    {
        MArray mx = mG.get_field("w");
        return mx.get_data<TWeight>();
    }
    
    const graph_size_t *o_degs() const
    {
        return mG.get_field("o_degs").get_data<graph_size_t>();
    }
    
    const graph_size_t *o_offsets() const
    {
        return mG.get_field("o_offsets").get_data<graph_size_t>();
    }
    
    const graph_size_t *i_degs() const
    {
        return mG.get_field("i_degs").get_data<graph_size_t>();
    }
    
    const graph_size_t *i_offsets() const
    {
        return mG.get_field("i_offsets").get_data<graph_size_t>();
    }
    

    // conversion to ref graph class
    
    inline CRefEdgeList<no_edge_weight> to_cref_edgelist() const
    {    
        return CRefEdgeList<no_edge_weight>(n(), m(), s(), t(), 0);
    } 
    
    template<typename TWeight>
    inline CRefEdgeList<TWeight> to_cref_wedgelist() const
    {
        return CRefEdgeList<TWeight>(n(), m(), s(), t(), w<TWeight>());
    }  
    
    inline CRefAdjList<no_edge_weight> to_cref_adjlist() const
    {
        return CRefAdjList<no_edge_weight>(n(), m(), s(), t(), 0, 
                o_degs(), o_offsets());
    }
    
    
    template<typename TWeight>
    inline CRefAdjList<TWeight> to_cref_wadjlist() const
    {
        return CRefAdjList<TWeight>(n(), m(), s(), t(), w<TWeight>(), 
                o_degs(), o_offsets());
    }
};


template<typename TWeight>
inline mxArray *create_matlab_edge_weight_matrix(int m, int n, const TWeight *src)
{
    return src_to_matlab_matrix<TWeight>(m, n, src);
}

template<>
inline mxArray *create_matlab_edge_weight_matrix<no_edge_weight>(
        int m, int n, const no_edge_weight*)
{
    return create_empty_matrix();
}




template<typename TWeight>
inline mxArray *create_matlab_graph_struct(
        graph_size_t n, graph_size_t m,
        const vertex_t *s, const vertex_t *t, const TWeight *w, 
        const graph_size_t *o_degs, const graph_size_t *o_offsets, 
        const graph_size_t *i_degs = 0, const graph_size_t *i_offsets = 0)
{
    static const char* fieldnames[] = {
        "n", "m", "s", "t", "w", 
        "o_degs", "o_offsets", "i_degs", "i_offsets"};
    static int nfields = 9;        
    
    mxArray *mxG = mxCreateStructMatrix(1, 1, nfields, fieldnames);
    
    mxSetField(mxG, 0, "n", create_matlab_scalar<double>(n));
    mxSetField(mxG, 0, "m", create_matlab_scalar<double>(m));
    
    mxSetField(mxG, 0, "s", src_to_matlab_matrix<graph_size_t>(m, 1, (const graph_size_t*)s));
    mxSetField(mxG, 0, "t", src_to_matlab_matrix<graph_size_t>(m, 1, (const graph_size_t*)t));        
    
    mxSetField(mxG, 0, "w", create_matlab_edge_weight_matrix<TWeight>(m, 1, w));
    
    mxSetField(mxG, 0, "o_degs", src_to_matlab_matrix_cond(n, 1, o_degs));
    mxSetField(mxG, 0, "o_offsets", src_to_matlab_matrix_cond(n, 1, o_offsets));
    mxSetField(mxG, 0, "i_degs", src_to_matlab_matrix_cond(n, 1, i_degs));
    mxSetField(mxG, 0, "i_offsets", src_to_matlab_matrix_cond(n, 1, i_offsets));
    
    return mxG;
}





}


#endif

