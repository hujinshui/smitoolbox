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

#include <string.h>
#include <vector>


using namespace smi;


char color_to_char(boost::default_color_type cr)
{
    if (cr == boost::color_traits<boost::default_color_type>::white())
    {
        return 'w';
    }
    else if (cr == boost::color_traits<boost::default_color_type>::gray())
    {
        return 'g';
    }
    else if (cr == boost::color_traits<boost::default_color_type>::black())
    {
        return 'b';
    }
    else
    {
        return 'o';
    }
}


struct vtoi
{
    typedef int result_type;
    
    int operator() (const vertex_t& v) const { return v.i + 1; }
};


/**
 * The struct to capture the information of BFS running
 */
struct BFSRecord
{
    typedef boost::default_color_type color_t;
    
    BFSRecord(graph_size_t n) 
    : colors(n, boost::color_traits<color_t>::white()), dist_map(0) { }
    
    ~BFSRecord()
    {
        if (dist_map != 0) delete[] dist_map;
    }
    
    
    std::vector<color_t> colors;
    
    std::vector<vertex_t> vertices;
    std::vector<vertex_t> parents;
    std::vector<int>      distances;  
    
    int *dist_map;
    
    
    void init_dist_map()
    {
        graph_size_t n = (graph_size_t)colors.size();
        dist_map = new int[n];
    }
    
    
    mxArray *vertices_to_matlab() const
    {
        return iter_to_matlab_row(vertices.begin(), vertices.size(), vtoi());
    }
    
    mxArray *parents_to_matlab() const
    {
        return iter_to_matlab_row(parents.begin(), parents.size(), vtoi());
    }
    
    mxArray *distances_to_matlab() const
    {
        return iter_to_matlab_row(distances.begin(), distances.size());
    }
    
    VertexRefMap<color_t> color_map() 
    {
        return &(colors[0]);
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
    
    void initialize_seed(const vertex_t& u) { }    
    
    void discover_vertex(const vertex_t& u, const CRefAdjList<no_edge_weight>& g)
    {
        m_record.vertices.push_back(u);
    }
    
private:
    BFSRecord& m_record;
};



/**
 * The BFS visitor that also records parents
 */
class BFSVisitorEx : public boost::default_bfs_visitor
{
public:    
    BFSVisitorEx(const CRefAdjList<no_edge_weight>& g, BFSRecord& rec) : m_record(rec) 
    {
        graph_size_t n = num_vertices(g);
        
        m_record.vertices.reserve(n);
        m_record.parents.reserve(n);        
    }
    
    void initialize_seed(const vertex_t& u)
    {
        m_record.vertices.push_back(u);
        m_record.parents.push_back(-1);
    }
        
    void tree_edge(const edge_t& e, const CRefAdjList<no_edge_weight>& g)
    {
        m_record.vertices.push_back(target(e, g));        
        m_record.parents.push_back(source(e, g));
    }        
    
private:
    BFSRecord& m_record;
};


/**
 * The BFS visitor that also records parents and distances
 */
class BFSVisitorExD : public boost::default_bfs_visitor
{
public:    
    BFSVisitorExD(const CRefAdjList<no_edge_weight>& g, BFSRecord& rec) : m_record(rec) 
    {
        graph_size_t n = num_vertices(g);
        
        m_record.vertices.reserve(n);
        m_record.parents.reserve(n); 
        m_record.distances.reserve(n);
        
        m_record.init_dist_map();
    }
    
    void initialize_seed(const vertex_t& u)
    {
        m_record.vertices.push_back(u);
        m_record.parents.push_back(-1);
        
        m_record.dist_map[u.i] = 0;
        m_record.distances.push_back(m_record.dist_map[u.i]);
    }
        
    void tree_edge(const edge_t& e, const CRefAdjList<no_edge_weight>& g)
    {
        vertex_t s = source(e, g);
        vertex_t t = target(e, g);
        
        m_record.vertices.push_back(t);        
        m_record.parents.push_back(s);
        
        m_record.dist_map[t.i] = m_record.dist_map[s.i] + 1;
        m_record.distances.push_back(m_record.dist_map[t.i]);
    }        
    
private:
    BFSRecord& m_record;    
};



/**
 * Core function
 */
template<class TVisitor>
void do_bfs(const CRefAdjList<no_edge_weight>& g, BFSRecord& record, int ns, int *s)
{
    using boost::visitor;
    using boost::color_map;
    
    typedef boost::default_color_type color_t;
    color_t white = boost::color_traits<color_t>::white();
    
    TVisitor vis(g, record);
    VertexRefMap<color_t> cmap = record.color_map();
    
        
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
            boost::breadth_first_visit(g, s[i], visitor(vis).color_map(cmap));
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
        do_bfs<BFSVisitorSimple>(g, record, ns, s);
        
        plhs[0] = record.vertices_to_matlab();
    }
    else if (nlhs == 2)
    {
        do_bfs<BFSVisitorEx>(g, record, ns, s);
        
        plhs[0] = record.vertices_to_matlab();
        plhs[1] = record.parents_to_matlab();
    }
    else if (nlhs == 3)
    {
        do_bfs<BFSVisitorExD>(g, record, ns, s);
        
        plhs[0] = record.vertices_to_matlab();
        plhs[1] = record.parents_to_matlab();
        plhs[2] = record.distances_to_matlab();
    }
    
}






