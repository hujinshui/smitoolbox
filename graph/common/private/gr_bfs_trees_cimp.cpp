/********************************************************************
 *
 *  gr_bfs_trees_cimp.cpp
 *
 *  The C++ mex implementation for gr_bfs.m
 *
 *  Created by Dahua Lin, on Oct 31, 2010
 *
 *******************************************************************/

#include "../../clib/graph_mex.h"

#include <boost/graph/breadth_first_search.hpp>

#include <vector>
#include <valarray>

using namespace smi;

/**
 * The struct to capture the information of BFS running
 */
struct BFSRecord
{
    typedef boost::default_color_type color_t;
    
    BFSRecord(graph_size_t n) 
    : num_vertices(n), colors(boost::white_color, n), parent_map(n), is_tree(true) { }
            
    // fields
    
    graph_size_t num_vertices;
    std::valarray<color_t> colors;
    std::valarray<vertex_t> parent_map;
    
    std::vector<vertex_t> vertices;
    bool is_tree;
        
    
    // actions
    
    void init_vertices()
    {
        vertices.reserve(num_vertices);
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
        
};



/**
 * The BFS visitor class
 */
class BFSTreeVisitor : public boost::default_bfs_visitor
{
public:    
    BFSTreeVisitor(BFSRecord& rec) 
    : vertices(rec.vertices), parent_map(rec.parent_map), is_tree(rec.is_tree)
    {           
    }
    
    void initialize_seed(const vertex_t& u)
    {
        vertices.push_back(u);
        parent_map[u.i] = -1;
    }
        
    void non_tree_edge(const edge_t& e, const CRefAdjList<no_edge_weight>& g)
    {
        vertex_t s = source(e, g);
        vertex_t t = target(e, g);
                
        if (t != parent_map[s.i]) is_tree = false;
    }    
    
    void tree_edge(const edge_t& e, const CRefAdjList<no_edge_weight>& g)
    {
        vertex_t s = source(e, g);
        vertex_t t = target(e, g);
        
        vertices.push_back(t);        
        parent_map[t.i] = s;
    }        
        
    
private:
    std::vector<vertex_t>& vertices;
    std::valarray<vertex_t>& parent_map;
    bool& is_tree;
};




/**
 * Core function
 */
void do_bfs_tree(const CRefAdjList<no_edge_weight>& g, BFSRecord& record, int ns, int *s)
{
    using boost::visitor;
    using boost::color_map;
    
    typedef boost::default_color_type color_t;
    color_t white = boost::color_traits<color_t>::white();
    
    BFSTreeVisitor vis(record);
    VertexRefMap<color_t> cmap = &(record.colors[0]);
                
    // subsequent search 
    
    boost::queue<vertex_t> Q;
    
    for (int i = 0; i < ns; ++i)
    {
        vertex_t v = s[i];
                        
        if (cmap[v] == white)
        {
            vis.initialize_seed(v);
            boost::breadth_first_visit(g, s[i], Q, vis, cmap);
        }
    }    
}



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
       
    record.init_vertices();
        
    do_bfs_tree(g, record, ns, s);
        
    plhs[0] = record.vertices_to_matlab();
    plhs[1] = record.parents_to_matlab();
    plhs[2] = mxCreateLogicalScalar(record.is_tree);
}






