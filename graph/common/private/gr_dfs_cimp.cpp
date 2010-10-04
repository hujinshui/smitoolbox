/********************************************************************
 *
 *  gr_dfs_cimp.cpp
 *
 *  The C++ mex implementation of gr_dfs
 *
 *  Created by Dahua Lin, on Oct 3, 2010
 *
 *******************************************************************/


#include "../../clib/mgraph.h"
#include "../../clib/graph_search.h"

using namespace smi; 


struct vnode_getter
{
    typedef int value_type;
    
    const SeqList<int>& m_vs;
    
    vnode_getter(const SeqList<int>& vs) : m_vs(vs) { }
    
    int count() const { return m_vs.size(); }
    
    int operator() (int i) const { return m_vs[i] + 1; }    
};


struct vpred_getter
{
    typedef int value_type;
    
    const SeqList<int>& m_vs;
    const DFSIteratorEx& m_it;
    
    vpred_getter(const SeqList<int>& vs, const DFSIteratorEx& bfs_it)
    : m_vs(vs), m_it(bfs_it) { }
    
    int count() const { return m_vs.size(); }
    
    int operator() (int i) const { return m_it.predecessor_of(m_vs[i]) + 1; }     
};


struct vdtime_getter
{
    typedef int value_type;
    
    const SeqList<int>& m_vs;
    const DFSIteratorEx& m_it;
    
    vdtime_getter(const SeqList<int>& vs, const DFSIteratorEx& bfs_it)
    : m_vs(vs), m_it(bfs_it) { }
    
    int count() const { return m_vs.size(); }
    
    int operator() (int i) const { return m_it.discover_time_of(m_vs[i]); }     
};


struct vftime_getter
{
    typedef int value_type;
    
    const SeqList<int>& m_vs;
    const DFSIteratorEx& m_it;
    
    vftime_getter(const SeqList<int>& vs, const DFSIteratorEx& bfs_it)
    : m_vs(vs), m_it(bfs_it) { }
    
    int count() const { return m_vs.size(); }
    
    int operator() (int i) const { return m_it.finish_time_of(m_vs[i]); }     
};



template<typename GIter>    
void do_dfs(GIter& dfs_it, int ns, const int *seeds, SeqList<int>& vs)
{    
    for (int i = 0; i < ns; ++i) 
    {
        int s = seeds[i];
        if (!dfs_it.is_discovered(s))
        {
            dfs_it.set_seed(s);
            int v;
            while ((v = dfs_it.next()) >= 0) 
            {
                vs.add(v);
            }
        }
    }   
}
       

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const MArray mG(prhs[0]);
    const MArray mSeeds(prhs[1]);
    
    // make graph
    
    RefGraph G = to_refgraph(mG);
    AdjList adjList(G);
    
    int ns = mSeeds.nelems();
    const int *seeds = mSeeds.get_data<int>();
    
    SeqList<int> vs(G.nnodes());
    
    // main
    
    if (nlhs <= 1)
    {
        DFSIterator dfs_it(adjList);
        do_dfs(dfs_it, ns, seeds, vs);
        plhs[0] = iter_to_matlab_row(vnode_getter(vs));
    }
    else if (nlhs <= 5)
    {
        DFSIteratorEx dfs_it(adjList);
        do_dfs(dfs_it, ns, seeds, vs);
        
        plhs[0] = iter_to_matlab_row(vnode_getter(vs));
        plhs[1] = iter_to_matlab_row(vpred_getter(vs, dfs_it));
        
        if (nlhs >= 3)
            plhs[2] = iter_to_matlab_row(vnode_getter(dfs_it.finish_order()));
        
        if (nlhs >= 4)
            plhs[3] = iter_to_matlab_row(vdtime_getter(vs, dfs_it));
        
        if (nlhs >= 5)
            plhs[4] = iter_to_matlab_row(vftime_getter(vs, dfs_it));        
                
    }

}
   



