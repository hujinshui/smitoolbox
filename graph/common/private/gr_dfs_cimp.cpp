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

#include <boost/graph/depth_first_search.hpp>

#include <vector>
#include <valarray>

using namespace bcs;
using namespace bcs::matlab;


/**
 * The struct to capture the information of DFS running
 */
struct DFSRecord
{
    typedef boost::default_color_type color_t;
    
    DFSRecord(gr_size_t n) 
    : num_vertices(n), colors(boost::white_color, n) 
    { }
        
    // fields
    
    gr_size_t num_vertices;
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
    
    marray vertices_to_matlab() const
    {
        size_t n = discover_seq.size();
        marray mA = create_marray<gr_index_t>(1, n);
        vertex_t *a = mA.data<vertex_t>();
        
        for (size_t i = 0; i < n; ++i)
        {
            a[i] = discover_seq[i].index + 1;
        }
        return mA;
    }
    
    marray parents_to_matlab() const
    {
        size_t n = parents_seq.size();
        marray mA = create_marray<gr_index_t>(1, n);
        vertex_t *a = mA.data<vertex_t>();
        
        for (size_t i = 0; i < n; ++i)
        {
            a[i] = parents_seq[i].index + 1;
        }
        return mA;
    }

    marray finishord_to_matlab() const
    {
        size_t n = finish_seq.size();
        marray mA = create_marray<gr_index_t>(1, n);
        vertex_t *a = mA.data<vertex_t>();
        
        for (size_t i = 0; i < n; ++i)
        {
            a[i] = finish_seq[i].index + 1;
        }
        return mA;       
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
    
    template<typename Graph>
    void discover_vertex(const vertex_t& u, const Graph& g)
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
        
    template<typename Graph>
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
        
    template<typename Graph>
    void tree_edge(const edge_t& e, const Graph& g)
    {
        vertex_t s = source(e, g);
        vertex_t t = target(e, g);
        
        vertices.push_back(t);        
        parents.push_back(s);         
    }        
    
    template<typename Graph>
    void finish_vertex(const vertex_t& u, const Graph& g)
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
template<typename TDir, class TVisitor>
void do_dfs(const_mgraph mG, const_aview1d<vertex_t> seeds, DFSRecord& record)
{
    gr_adjlist<TDir> g = to_gr_adjlist<TDir>(mG); 
        
    typedef boost::default_color_type color_t;
    
    TVisitor vis(record);
    vertex_ref_map<color_t> cmap( &(record.colors[0]) );
            
    // initial search with s[0]
    
    for (index_t i = 0; i < (index_t)seeds.nelems(); ++i)
    {
        const vertex_t& v = seeds[i];
                        
        if (cmap[v] == boost::white_color)
        {
            vis.initialize_seed(v);
            boost::depth_first_visit(g, v, vis, cmap);
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
    // take input
    
    const_mgraph mG(prhs[0]);
    const_aview1d<vertex_t> seeds = view1d<vertex_t>(prhs[1]);
    
    char dtype = mG.dtype();
    
    // main
    
    DFSRecord record(mG.nv());    
       
    if (nlhs <= 1)
    {        
        record.init_vertices();
        
        if (dtype == 'd')
        {
            do_dfs<gr_directed, DFSVisitorSimple>(mG, seeds, record);            
        }
        else if (dtype == 'u')
        {
            do_dfs<gr_undirected, DFSVisitorSimple>(mG, seeds, record);  
        }
    }
    else if (nlhs == 2)
    {
        record.init_vertices();
        record.init_parents();
        
        if (dtype == 'd')
        {
            do_dfs<gr_directed, DFSVisitorEx>(mG, seeds, record);
        }
        else if (dtype == 'u')
        {
            do_dfs<gr_undirected, DFSVisitorEx>(mG, seeds, record);  
        }
    }
    else if (nlhs == 3)
    {
        record.init_vertices();
        record.init_parents();
        record.init_finish_order();
        
        if (dtype == 'd')
        {
            do_dfs<gr_directed, DFSVisitorExF>(mG, seeds, record);
        }
        else if (dtype == 'u')
        {
            do_dfs<gr_undirected, DFSVisitorExF>(mG, seeds, record);
        }
    }
        
    plhs[0] = record.vertices_to_matlab().mx_ptr();
    if (nlhs >= 2)
        plhs[1] = record.parents_to_matlab().mx_ptr();
    if (nlhs >= 3)
        plhs[2] = record.finishord_to_matlab().mx_ptr();    
}






