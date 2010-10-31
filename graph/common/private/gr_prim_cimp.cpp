/********************************************************************
 *
 *  gr_prim_cimp.cpp
 *
 *  The C++ mex implementation for Prim's minimum spanning tree
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 ********************************************************************/


#include "../../clib/graph_mex.h"

#include <vector>
#include <valarray>
#include <iterator>

#include <boost/graph/prim_minimum_spanning_tree.hpp>


typedef boost::default_color_type color_t;


using namespace smi;

template<typename TWeight>
struct PrimRecord
{
    std::vector<vertex_t> vertices;     // the vertices in finished order
    
    std::valarray<vertex_t> parent_map;   
    std::valarray<TWeight> eweight_map;
    std::valarray<color_t> color_map;
    std::valarray<graph_size_t> edge_map;
    
    PrimRecord(graph_size_t n)
    : parent_map(n)
    , eweight_map(n)
    , color_map(boost::color_traits<color_t>::white(), n)
    , edge_map(n)
    {
        vertices.reserve(n);
    }        
       
    mxArray *vertices_to_matlab() const
    {
        return iter_to_matlab_column(++vertices.begin(), vertices.size()-1, 
                vertex_to_mindex());        
    }
    
    mxArray *edges_to_matlab() const
    {
        return iter_to_matlab_column(++vertices.begin(), vertices.size()-1,
                unary_chain(vertex_to_index(), arr_map(edge_map), vertex_to_mindex()) );                
    } 
    
    mxArray *parents_to_matlab() const
    {
        return iter_to_matlab_column(++vertices.begin(), vertices.size()-1,
                unary_chain(vertex_to_index(), arr_map(parent_map), vertex_to_mindex()) );                
    }    
    
    mxArray *eweights_to_matlab() const
    {
        return iter_to_matlab_column(++vertices.begin(), vertices.size()-1,
                unary_chain(vertex_to_index(), arr_map(eweight_map)) );                
    }         
    
};


template<typename TWeight>
class PrimVisitor : public boost::default_dijkstra_visitor
{
public:
    
    PrimVisitor(PrimRecord<TWeight>& record)
    : m_record(record)
    {        
    }
    
    void edge_relaxed(const edge_t& e, const CRefAdjList<TWeight, boost::undirected_tag>& g)
    {
        vertex_t s = source(e, g);
        vertex_t t = target(e, g);
        
        m_record.edge_map[t.i] = e.i;
    }
    
    void finish_vertex(const vertex_t& u, const CRefAdjList<TWeight, boost::undirected_tag>& g)
    {
        m_record.vertices.push_back(u);                
    }    
    
private:
    PrimRecord<TWeight>& m_record;
};


template<typename TWeight>
void main_delegate(const matlab_graph_repr& gr, int s, int nlhs, mxArray *plhs[])
{
    CRefAdjList<TWeight, boost::undirected_tag> g = gr.to_cref_wadjlist_ud<TWeight>();
    
    PrimRecord<TWeight> record(num_vertices(g));    
    PrimVisitor<TWeight> vis(record);

    using boost::root_vertex;
    using boost::distance_map;
    using boost::color_map;
    
    VertexRefMap<vertex_t> pmap = &(record.parent_map[0]);
    VertexRefMap<TWeight> ewmap = &(record.eweight_map[0]);
    VertexRefMap<color_t> cmap = &(record.color_map[0]);
        
    vertex_t rv(s);
    
    boost::prim_minimum_spanning_tree(g, pmap, 
             root_vertex(rv).distance_map(ewmap).color_map(cmap).visitor(vis));
    
    pmap[0] = -1;
    
    plhs[0] = record.edges_to_matlab();
}



/***
 * main entry
 *
 * Inputs:
 *  [0]: G:     the undirected graph with adjlist representation
 *  [1]: r:     the root node of the MST
 *
 * Outputs:
 *  [0]: s:     the sources of edges
 *  [1]: t:     the targets of edges
 *  [2]: w:     the weights of edges
 *  [3]: d:     the distances from each end-node to root
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    matlab_graph_repr gr(prhs[0]);
    int s = MArray(prhs[1]).get_scalar<int>();
    
    switch (gr.weight_class())
    {
        case mxDOUBLE_CLASS:
            main_delegate<double>(gr, s, nlhs, plhs);
            break;
            
        case mxSINGLE_CLASS:
            main_delegate<float>(gr, s, nlhs, plhs);
            break;
            
        case mxINT32_CLASS:
            main_delegate<int>(gr, s, nlhs, plhs);
            break;
            
        case mxUINT32_CLASS:
            main_delegate<unsigned int>(gr, s, nlhs, plhs);
            break;
            
        default:
            mexErrMsgIdAndTxt("gr_prim_mst:invalidarg", 
                    "The weight value should be double, single, int32, or uint32.");
    }
}


