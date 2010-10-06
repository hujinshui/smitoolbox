/********************************************************************
 *
 *  gr_bfs_cimp.cpp
 *
 *  The C++ mex implementation of gr_bfs
 *
 *  Created by Dahua Lin, on Oct 3, 2010
 *
 *******************************************************************/


#include "../../clib/mgraph.h"
#include "../../clib/graph_search.h"

using namespace smi; 


// Getters 

struct node_getter
{    
    int operator() (int v) const 
    { 
        return v + 1; 
    }    
};

template<typename TAlgIter>
struct PredGetter
{
    const TAlgIter& m_it;    
    PredGetter(const TAlgIter& it) : m_it(it) { }
        
    int operator() (int v) const 
    { 
        return m_it.predecessor_of(v) + 1; 
    }     
};

template<typename TAlgIter>
PredGetter<TAlgIter> pred_getter(const TAlgIter& it)
{
    return PredGetter<TAlgIter>(it);
}

template<typename TAlgIter>
struct DistGetter
{
    const TAlgIter& m_it;    
    DistGetter(const TAlgIter& it) : m_it(it) { }
        
    int operator() (int v) const 
    { 
        return m_it.distance_of(v); 
    } 
};

template<typename TAlgIter>
DistGetter<TAlgIter> dist_getter(const TAlgIter& it)
{
    return DistGetter<TAlgIter>(it);
}

template<typename TAlgIter>
struct DTimeGetter
{
    const TAlgIter& m_it;    
    DTimeGetter(const TAlgIter& it) : m_it(it) { }
        
    int operator() (int v) const 
    { 
        return m_it.discover_time_of(v); 
    } 
};

template<typename TAlgIter>
DTimeGetter<TAlgIter> dtime_getter(const TAlgIter& it)
{
    return DTimeGetter<TAlgIter>(it);
}

template<typename TAlgIter>
struct FTimeGetter
{
    const TAlgIter& m_it;    
    FTimeGetter(const TAlgIter& it) : m_it(it) { }
        
    int operator() (int v) const 
    { 
        return m_it.finish_time_of(v); 
    } 
};

template<typename TAlgIter>
FTimeGetter<TAlgIter> ftime_getter(const TAlgIter& it)
{
    return FTimeGetter<TAlgIter>(it);
}





// Core algorithm

template<typename TAlgIter>    
void do_bfs(TAlgIter& bfs_it, int ns, const int *seeds, std::vector<int>& vs)
{    
    for (int i = 0; i < ns; ++i) 
        bfs_it.add_seed(seeds[i]);
        
    int v = -1;
    while ((v = bfs_it.next()) >= 0)
    {
        vs.push_back(v);
    }    
    
    int n = vs.size();
}
       

template<typename TAlgIter>    
void do_dfs(TAlgIter& dfs_it, int ns, const int *seeds, std::vector<int>& vs)
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
                vs.push_back(v);
            }
        }
    }   
}





// Main Entry

    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const MArray mG(prhs[0]);
    const MArray mSeeds(prhs[1]);
    const MArray mCode(prhs[2]);
    
    char code = (char)mCode.get_scalar<mxChar>();    
    
    // make graph
    
    RefGraph G = to_refgraph(mG);
    AdjList adjList(G);
    
    int ns = mSeeds.nelems();
    const int *seeds = mSeeds.get_data<int>();
    
    std::vector<int> vs;
    vs.reserve(G.nnodes());
    
    // main
    
    if (code == 'b')
    {
        if (nlhs <= 1)
        {
            BFSIterator bfs_it(adjList);
            do_bfs(bfs_it, ns, seeds, vs);        

            plhs[0] = iter_to_matlab_row(vs.begin(), vs.size(), node_getter());
        }    
        else
        {
            BFSIteratorEx bfs_it(adjList);
            do_bfs(bfs_it, ns, seeds, vs);

            plhs[0] = iter_to_matlab_row(vs.begin(), (int)vs.size(), node_getter());
            plhs[1] = iter_to_matlab_row(vs.begin(), (int)vs.size(), pred_getter(bfs_it));
            if (nlhs >= 3)
                plhs[2] = iter_to_matlab_row(vs.begin(), (int)vs.size(), dist_getter(bfs_it));
        }
    }
    else if (code == 'd')
    {
        if (nlhs <= 1)
        {
            DFSIterator dfs_it(adjList);
            do_dfs(dfs_it, ns, seeds, vs);
            plhs[0] = iter_to_matlab_row(vs.begin(), (int)vs.size(), node_getter());
        }
        else if (nlhs <= 5)
        {
            DFSIteratorEx dfs_it(adjList);
            do_dfs(dfs_it, ns, seeds, vs);

            plhs[0] = iter_to_matlab_row(vs.begin(), (int)vs.size(), node_getter());
            plhs[1] = iter_to_matlab_row(vs.begin(), (int)vs.size(), pred_getter(dfs_it));

            if (nlhs >= 3)
            {
                const std::vector<int>& ford = dfs_it.finish_order();
                plhs[2] = iter_to_matlab_row(ford.begin(), (int)ford.size(), node_getter());
            }

            if (nlhs >= 4)
                plhs[3] = iter_to_matlab_row(vs.begin(), (int)vs.size(), dtime_getter(dfs_it));

            if (nlhs >= 5)
                plhs[4] = iter_to_matlab_row(vs.begin(), (int)vs.size(), ftime_getter(dfs_it)); 
        }
    }
}
   



