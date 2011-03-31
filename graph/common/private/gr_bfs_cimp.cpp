/********************************************************************
 *
 *  gr_bfs_cimp.cpp
 *
 *  The C++ mex implementation for gr_bfs.m
 *
 *  Created by Dahua Lin, on Oct 9, 2010
 *
 *******************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/matlab/mgraph.h>
#include <bcslib/graph/bgl_port.h>

#include <boost/graph/breadth_first_search.hpp>

#include <vector>
#include <valarray>

using namespace bcs;
using namespace bcs::matlab;

/**
 * The struct to capture the information of BFS running
 */
struct BFSRecord
{
    typedef boost::default_color_type color_t;
    
    // fields
    
    gr_size_t num_vertices;
    std::valarray<color_t> colors;
    
    std::vector<vertex_t> vertices;
    std::vector<vertex_t> parents;    
    std::valarray<int32_t> *p_dist_map;
    
    // constructors
    
    BFSRecord(gr_size_t n) 
    : num_vertices(n), colors(boost::white_color, n), p_dist_map(0) { }
    
    ~BFSRecord()
    {
        if (p_dist_map != 0) delete p_dist_map;
    }                
    
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
    
    
    marray vertices_to_matlab() const
    {
        size_t n = vertices.size();
        marray mA = create_marray<gr_index_t>(1, n);
        vertex_t *a = mA.data<vertex_t>();
        
        for (size_t i = 0; i < n; ++i)
        {
            a[i] = vertices[i].index + 1;
        }
        return mA;
    }
    
    marray parents_to_matlab() const
    {
        size_t n = parents.size();
        marray mA = create_marray<gr_index_t>(1, n);
        vertex_t *a = mA.data<vertex_t>();
        
        for (size_t i = 0; i < n; ++i)
        {
            a[i] = parents[i].index + 1;
        }
        return mA;
    }
    
    marray distances_to_matlab() const
    {
        const std::valarray<int>& dist_map = *p_dist_map;
        
        size_t n = vertices.size();
        marray mA = create_marray<int32_t>(1, n);
        int32_t *a = mA.data<int32_t>();
        
        for (size_t i = 0; i < n; ++i)
        {
            a[i] = dist_map[vertices[i].index];
        }
        return mA;
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
    
    template<class Graph>
    void discover_vertex(const vertex_t& u, const Graph& g)
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
        
    template<class Graph>
    void tree_edge(const edge_t& e, const Graph& g)
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
        dists[u.index] = 0;
    }
        
    template<class Graph>
    void tree_edge(const edge_t& e, const Graph& g)
    {
        vertex_t s = source(e, g);
        vertex_t t = target(e, g);
        
        vertices.push_back(t);        
        parents.push_back(s);        
        dists[t.index] = dists[s.index] + 1;
    }        
    
private:
    std::vector<vertex_t>& vertices;
    std::vector<vertex_t>& parents;
    std::valarray<int>& dists;
};



// core function
template<typename TDir, class TVisitor>
void do_bfs(const_mgraph mG, const_aview1d<vertex_t> seeds, BFSRecord& record)
{
    // take graph
    
    gr_adjlist<TDir> g = to_gr_adjlist<TDir>(mG);
    
    // prepare params
    
    typedef boost::default_color_type color_t;
    
    boost::queue<vertex_t> Q;
    TVisitor vis(record);
    vertex_ref_map<color_t> cmap( &(record.colors[0]) );
    
    for (index_t i = 0; i < (index_t)seeds.nelems(); ++i)
    {
        const vertex_t& v = seeds[i];
        
        if (cmap[v] == boost::white_color)
        {
            vis.initialize_seed(v);
            boost::breadth_first_visit(g, v, Q, vis, cmap);
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
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    // take input
    
    const_mgraph mG(prhs[0]);
    const_aview1d<vertex_t> seeds = view1d<vertex_t>(prhs[1]);
    
    char dtype = mG.dtype();
    
    // main
    BFSRecord record(mG.nv());    
       
    if (nlhs <= 1)
    {        
        record.init_vertices();
        
        if (dtype == 'd')
        {
            do_bfs<gr_directed, BFSVisitorSimple>(mG, seeds, record);
        }
        else if (dtype == 'u')
        {
            do_bfs<gr_undirected, BFSVisitorSimple>(mG, seeds, record);
        }
    }
    else if (nlhs == 2)
    {
        record.init_vertices();
        record.init_parents();
        
        if (dtype == 'd')
        {
            do_bfs<gr_directed, BFSVisitorEx>(mG, seeds, record);
        }
        else if (dtype == 'u')
        {
            do_bfs<gr_undirected, BFSVisitorEx>(mG, seeds, record);
        }
    }
    else if (nlhs == 3)
    {
        record.init_vertices();
        record.init_parents();
        record.init_distmap();
        
        if (dtype == 'd')
        {
            do_bfs<gr_directed, BFSVisitorExD>(mG, seeds, record); 
        }
        else if (dtype == 'u')
        {
            do_bfs<gr_undirected, BFSVisitorExD>(mG, seeds, record); 
        }
    }
        
    plhs[0] = record.vertices_to_matlab().mx_ptr();
    if (nlhs >= 2)
        plhs[1] = record.parents_to_matlab().mx_ptr();
    if (nlhs >= 3)
        plhs[2] = record.distances_to_matlab().mx_ptr();    
}


BCSMEX_MAINDEF



