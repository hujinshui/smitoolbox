/********************************************************************
 *
 *  gr_bfs_cimp.cpp
 *
 *  The C++ mex implementation for gr_bfs.m
 *
 *  Created by Dahua Lin, on Oct 9, 2010
 *
 *******************************************************************/

#include "../../clib/graph_mex.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <vector>


using namespace smi;


struct vtoi
{
    typedef int result_type;
    
    int operator() (const vertex_t& v) const { return v.i; }
};







/**
 * The struct to capture the information of BFS running
 */
struct BFSRecord
{
    BFSRecord(graph_size_t n) : colors(n) { }
    
    VertexColorMap colors;
    
    std::vector<vertex_t> vertices;
    std::vector<vertex_t> parents;
    std::vector<int>      distances;
    
    mxArray *vertices_to_matlab() const
    {
        return iter_to_matlab_row(vertices.begin(), vertices.size(), vtoi());
    }
    
};



/**
 * The BFS visitor that only records discovered vertices
 */
class BFSVisitorSimple : public boost::default_bfs_visitor
{
public:    
    BFSVisitorSimple(const CRefAdjList<no_edge_weight>& g, BFSRecord& rec) : m_record(rec) 
    {
        m_record.vertices.reserve(num_vertices(g));
    }
    
    void discover_vertex(const vertex_t& u, const CRefAdjList<no_edge_weight>& g)
    {
        m_record.vertices.push_back(u);
    }
    
private:
    BFSRecord& m_record;
};


/***
 * main entry
 *
 * Inputs:
 *  [0]: G:     the graph with adjlist representation
 *  [1]: s:     the source node(s)
 *
 * Outputs:
 *  [0]: vs:        the vertices in discovery order
 *  [1]: parents:   parents in the search tree (corresponding to vs)
 *  [2]: dists:     the distances to the root in search tree (corresponding to vs)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CRefAdjList<no_edge_weight> g = matlab_graph_repr(prhs[0]).to_cref_adjlist();
    
    MArray mS(prhs[1]);
    int ns = mS.nelems();
    int *s = mS.get_data<int>();
    
    BFSRecord record(num_vertices(g));    
    
    using boost::visitor;
    using boost::color_map;
    
    boost::adjacency_list<>* padj = 0;
    boost::adjacency_list<>::vertex_descriptor v0;
    
    if (nlhs <= 1)
    {
        BFSVisitorSimple vis(g, record);
        
        vertex_t v0(s[0]);
        boost::breadth_first_search(g, v0, 
                visitor(boost::default_bfs_visitor()).color_map(record.colors));                
        
        plhs[0] = record.vertices_to_matlab();
    }    
    
}






