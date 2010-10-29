/********************************************************************
 *
 *  gr_dijkstra_cimp.cpp
 *
 *  The C++ mex implementation for Dijkstra shortest-path solving
 *  algorithm
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 ********************************************************************/

#include "../../clib/graph_mex.h"

#include <vector>
#include <valarray>

#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace smi;

typedef boost::default_color_type color_t;


template<typename TWeight>
struct DijkstraRecord
{
    std::vector<vertex_t> vertices;     // the vertices in finished order
    
    std::valarray<vertex_t> parent_map;                        
    std::valarray<TWeight> distance_map;
    std::valarray<color_t> color_map;
    
    DijkstraRecord(graph_size_t n)
    : parent_map(n)
    , distance_map(TWeight(0), n)
    , color_map(boost::color_traits<color_t>::white(), n)
    {
        vertices.reserve(n);
    }        
       
    mxArray *vertices_to_matlab() const
    {
        return iter_to_matlab_row(vertices.begin(), vertices.size(), 
                vertex_to_mindex());
        
    }
    
    mxArray *parents_to_matlab() const
    {
        return iter_to_matlab_row(vertices.begin(), vertices.size(),
                unary_chain(vertex_to_index(), arr_map(parent_map), vertex_to_mindex()) );                
    }
    
    
    mxArray *distances_to_matlab() const
    {
        return iter_to_matlab_row(vertices.begin(), vertices.size(),
                unary_chain(vertex_to_index(), arr_map(distance_map)) );                
    }         
    
};


template<typename TWeight>
class DijkstraVisitorSimple : public boost::default_dijkstra_visitor
{
public:
    
    DijkstraVisitorSimple(DijkstraRecord<TWeight>& record)
    : m_record(record)
    {        
    }
    
    void finish_vertex(const vertex_t& u, const CRefAdjList<TWeight>& g)
    {
        m_record.vertices.push_back(u);
    }    
    
private:
    DijkstraRecord<TWeight>& m_record;
};





template<typename TWeight>
void main_delegate(const matlab_graph_repr& gr, int s, int nlhs, mxArray *plhs[])
{
    CRefAdjList<TWeight> g = gr.to_cref_wadjlist<TWeight>();
    
    DijkstraRecord<TWeight> record(num_vertices(g));    
    DijkstraVisitorSimple<TWeight> vis(record);
    
    using boost::weight_map;
    using boost::predecessor_map;
    using boost::distance_map;
    using boost::color_map;
    
    VertexRefMap<vertex_t> pmap = &(record.parent_map[0]);
    VertexRefMap<TWeight> dmap = &(record.distance_map[0]);
    VertexRefMap<color_t> cmap = &(record.color_map[0]);
        
    boost::dijkstra_shortest_paths(g, vertex_t(s), 
            predecessor_map(pmap).distance_map(dmap).color_map(cmap).visitor(vis));
    
    pmap[0] = -1;
    
    plhs[0] = record.vertices_to_matlab();
    if (nlhs >= 2)
        plhs[1] = record.parents_to_matlab();
    if (nlhs >= 3)
        plhs[2] = record.distances_to_matlab();
}



/***
 * main entry
 *
 * Inputs:
 *  [0]: G:     the graph with adjlist representation
 *  [1]: s:     the source node
 *
 * Outputs:
 *  [0]: vs:        the vertices in discovery order
 *  [1]: parents:   parents in the search tree (corresponding to vs)
 *  [2]: dists:     the distances to the source in search tree (corresponding to vs)
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
            mexErrMsgIdAndTxt("gr_dijkstra_cimp:invalidarg", 
                    "The weight value should be double, single, int32, or uint32.");
    }
}
