/********************************************************************
 *
 *  gr_bfs_cimp.cpp
 *
 *  The C++ mex implementation of gr_bfs_cimp
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
    const BFSIteratorEx& m_it;
    
    vpred_getter(const SeqList<int>& vs, const BFSIteratorEx& bfs_it)
    : m_vs(vs), m_it(bfs_it) { }
    
    int count() const { return m_vs.size(); }
    
    int operator() (int i) const { return m_it.predecessor_of(m_vs[i]) + 1; }     
};


struct vdist_getter
{
    typedef int value_type;
    
    const SeqList<int>& m_vs;
    const BFSIteratorEx& m_it;
    
    vdist_getter(const SeqList<int>& vs, const BFSIteratorEx& bfs_it)
    : m_vs(vs), m_it(bfs_it) { }
    
    int count() const { return m_vs.size(); }
    
    int operator() (int i) const { return m_it.distance_of(m_vs[i]); }     
};



template<typename GIter>    
void do_bfs(GIter& bfs_it, int ns, const int *seeds, SeqList<int>& vs)
{    
    for (int i = 0; i < ns; ++i) 
        bfs_it.add_seed(seeds[i]);
        
    int v = -1;
    while ((v = bfs_it.next()) >= 0)
    {
        vs.add(v);
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
        BFSIterator bfs_it(adjList);
        do_bfs(bfs_it, ns, seeds, vs);
        plhs[0] = iter_to_matlab_row(vnode_getter(vs));
    }    
    else
    {
        BFSIteratorEx bfs_it(adjList);
        do_bfs(bfs_it, ns, seeds, vs);
        
        plhs[0] = iter_to_matlab_row(vnode_getter(vs));
        plhs[1] = iter_to_matlab_row(vpred_getter(vs, bfs_it));
        if (nlhs >= 3)
            plhs[2] = iter_to_matlab_row(vdist_getter(vs, bfs_it));
    }
}
   



