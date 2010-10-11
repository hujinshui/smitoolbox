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

#include <boost/graph/depth_first_search.hpp>

#include <vector>
#include <valarray>


using namespace smi;

/**
 * The struct to capture the information of DFS running
 */
struct DFSRecord
{
    typedef boost::default_color_type color_t;
    
    DFSRecord(graph_size_t n) 
    : num_vertices(n), colors(boost::white_color, n) 
    { }
        
    // fields
    
    graph_size_t num_vertices;
    std::valarray<color_t> colors;
    
    std::vector<vertex_t> discover_seq;
    std::vector<vertex_t> parents_seq;
    std::vector<vertex_t> finish_seq;
    
        
    // actions
    
    void init_vertices()
    {
        discover_seq.reserve(num_vertices);
    }
    
    void init_parents()
    {
        parents_seq.reserve(num_vertices);
    }
    
    void init_finish_order()
    {
        finish_seq.reserve(num_vertices);
    }
        
    // output
    
    mxArray *vertices_to_matlab() const
    {
        return iter_to_matlab_row(discover_seq.begin(), discover_seq.size(), 
                vertex_to_mindex());
    }
    
    mxArray *parents_to_matlab() const
    {
        return iter_to_matlab_row(parents_seq.begin(), parents_seq.size(), 
                vertex_to_mindex());
    }

    mxArray *finishord_to_matlab() const
    {
        return iter_to_matlab_row(finish_seq.begin(), finish_seq.size(),
                vertex_to_mindex());
    }
    
};


/**
 * The BFS visitor that only records discovered vertices
 */
class DFSVisitorSimple : public boost::default_dfs_visitor
{
public:    
    DFSVisitorSimple(DFSRecord& rec) 
    : vertices(rec.discover_seq) 
    {
    }
    
    void initialize_seed(const vertex_t& u) { }    
    
    void discover_vertex(const vertex_t& u, const CRefAdjList<no_edge_weight>& g)
    {
        vertices.push_back(u);
    }
    
private:
    std::vector<vertex_t>& vertices;
};


/**
 * The DFS visitor that also records parents
 */
class DFSVisitorEx : public boost::default_dfs_visitor
{
public:    
    DFSVisitorEx(DFSRecord& rec) 
    : vertices(rec.discover_seq), parents(rec.parents_seq)
    {           
    }
    
    void initialize_seed(const vertex_t& u)
    {
        vertices.push_back(u);
        parents.push_back(-1);
    }
        
    void tree_edge(const edge_t& e, const CRefAdjList<no_edge_weight>& g)
    {
        vertices.push_back(target(e, g));        
        parents.push_back(source(e, g));
    }        
    
private:
    std::vector<vertex_t>& vertices;
    std::vector<vertex_t>& parents;
};


/**
 * The DFS visitor that also records parents and finishing order
 */
class DFSVisitorExF : public boost::default_dfs_visitor
{
public:    
    DFSVisitorExF(DFSRecord& rec) 
    : vertices(rec.discover_seq), parents(rec.parents_seq), ford(rec.finish_seq)
    {        
    }
    
    void initialize_seed(const vertex_t& u)
    {
        vertices.push_back(u);
        parents.push_back(-1);         
    }
        
    void tree_edge(const edge_t& e, const CRefAdjList<no_edge_weight>& g)
    {
        vertex_t s = source(e, g);
        vertex_t t = target(e, g);
        
        vertices.push_back(t);        
        parents.push_back(s);         
    }        
    
    void finish_vertex(const vertex_t& u, const CRefAdjList<no_edge_weight>& g)
    {
        ford.push_back(u);
    }
    
private:
    std::vector<vertex_t>& vertices;
    std::vector<vertex_t>& parents;
    std::vector<vertex_t>& ford;
};



/**
 * Core function
 */
template<typename TVisitor>
void do_dfs(const CRefAdjList<no_edge_weight>& g, DFSRecord& record, int ns, int *s)
{
    using boost::visitor;
    using boost::color_map;
    
    typedef boost::default_color_type color_t;
    color_t white = boost::color_traits<color_t>::white();
    
    TVisitor vis(record);
    VertexRefMap<color_t> cmap = &(record.colors[0]);
            
    // initial search with s[0]
    
    for (int i = 0; i < ns; ++i)
    {
        vertex_t v = s[i];
                        
        if (cmap[v] == white)
        {
            vis.initialize_seed(v);
            boost::depth_first_visit(g, s[i], vis, cmap);
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
    
    DFSRecord record(num_vertices(g));    
       
    if (nlhs <= 1)
    {        
        record.init_vertices();
        
        do_dfs<DFSVisitorSimple>(g, record, ns, s);
    }
    else if (nlhs == 2)
    {
        record.init_vertices();
        record.init_parents();
        
        do_dfs<DFSVisitorEx>(g, record, ns, s);
    }
    else if (nlhs == 3)
    {
        record.init_vertices();
        record.init_parents();
        record.init_finish_order();
        
        do_dfs<DFSVisitorExF>(g, record, ns, s);                
    }
        
    plhs[0] = record.vertices_to_matlab();
    if (nlhs >= 2)
        plhs[1] = record.parents_to_matlab();
    if (nlhs >= 3)
        plhs[2] = record.finishord_to_matlab();    
}






