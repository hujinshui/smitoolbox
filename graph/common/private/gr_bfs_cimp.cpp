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
    : num_vertices(n), colors(boost::white_color, n), p_dist_map(0) { }
    
    ~BFSRecord()
    {
        if (p_dist_map != 0) delete p_dist_map;
    }
        
    // fields
    
    graph_size_t num_vertices;
    std::valarray<color_t> colors;
    
    std::vector<vertex_t> vertices;
    std::vector<vertex_t> parents;    
    std::valarray<int> *p_dist_map;
        
    
    // actions
    
    void init_vertices()
    {
        vertices.reserve(num_vertices);
    }
    
    void init_parents()
    {
        parents.reserve(num_vertices);
    }
    
    void init_distmap()
    {
        p_dist_map = new std::valarray<int>(0, num_vertices);
    }
    
    
    mxArray *vertices_to_matlab() const
    {
        return iter_to_matlab_row(vertices.begin(), vertices.size(), 
                vertex_to_mindex());
    }
    
    mxArray *parents_to_matlab() const
    {
        return iter_to_matlab_row(parents.begin(), parents.size(), 
                vertex_to_mindex());
    }
    
    mxArray *distances_to_matlab() const
    {
        return iter_to_matlab_row(vertices.begin(), vertices.size(), 
                unary_chain(vertex_to_index(), arr_map(*p_dist_map)));
    }   
    
};


/**
 * The BFS visitor that only records discovered vertices
 */
class BFSVisitorSimple : public boost::default_bfs_visitor
{
public:    
    BFSVisitorSimple(BFSRecord& rec) 
    : vertices(rec.vertices) 
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
 * The BFS visitor that also records parents
 */
class BFSVisitorEx : public boost::default_bfs_visitor
{
public:    
    BFSVisitorEx(BFSRecord& rec) 
    : vertices(rec.vertices), parents(rec.parents)
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
 * The BFS visitor that also records parents and distances
 */
class BFSVisitorExD : public boost::default_bfs_visitor
{
public:    
    BFSVisitorExD(BFSRecord& rec) 
    : vertices(rec.vertices), parents(rec.parents), dists(*(rec.p_dist_map))
    {        
    }
    
    void initialize_seed(const vertex_t& u)
    {
        vertices.push_back(u);
        parents.push_back(-1);        
        dists[u.i] = 0;
    }
        
    void tree_edge(const edge_t& e, const CRefAdjList<no_edge_weight>& g)
    {
        vertex_t s = source(e, g);
        vertex_t t = target(e, g);
        
        vertices.push_back(t);        
        parents.push_back(s);        
        dists[t.i] = dists[s.i] + 1;
    }        
    
private:
    std::vector<vertex_t>& vertices;
    std::vector<vertex_t>& parents;
    std::valarray<int>& dists;
};



/**
 * Core function
 */
template<typename TVisitor>
void do_bfs(const CRefAdjList<no_edge_weight>& g, BFSRecord& record, int ns, int *s)
{
    using boost::visitor;
    using boost::color_map;
    
    typedef boost::default_color_type color_t;
    color_t white = boost::color_traits<color_t>::white();
    
    TVisitor vis(record);
    VertexRefMap<color_t> cmap = &(record.colors[0]);
            
    // initial search with s[0]
    
    vis.initialize_seed(s[0]);
    boost::breadth_first_search(g, s[0], visitor(vis).color_map(cmap));
    
    // subsequent search 
    
    for (int i = 1; i < ns; ++i)
    {
        vertex_t v = s[i];
                        
        if (cmap[v] == white)
        {
            vis.initialize_seed(v);
            boost::breadth_first_visit(g, s[i], 
                    visitor(vis).color_map(cmap));
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
       
    if (nlhs <= 1)
    {        
        record.init_vertices();
        
        do_bfs<BFSVisitorSimple>(g, record, ns, s);
    }
    else if (nlhs == 2)
    {
        record.init_vertices();
        record.init_parents();
        
        do_bfs<BFSVisitorEx>(g, record, ns, s);
    }
    else if (nlhs == 3)
    {
        record.init_vertices();
        record.init_parents();
        record.init_distmap();
        
        do_bfs<BFSVisitorExD>(g, record, ns, s);                
    }
        
    plhs[0] = record.vertices_to_matlab();
    if (nlhs >= 2)
        plhs[1] = record.parents_to_matlab();
    if (nlhs >= 3)
        plhs[2] = record.distances_to_matlab();    
}






