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

#include <functional>
#include <algorithm>


// Auxiliary functions

template<typename T>
inline mxArray* v2mrow(const std::vector<T>& vec)        
{
    return iter_to_matlab_row(vec.begin(), (int)vec.size());
}



// Core algorithm

/***********************************************************
 *
 *  BFS 
 *
 **********************************************************/
      

struct bfs_observer2
{
    bfs_observer2(int n) : preds(n)
    {
    }
    
    void on_visit(int p, int v)
    {
        preds[v] = p;
    }
    
    void on_discover(int p, int v) { }
    
    Array<int> preds;
};


struct bfs_observer3
{
    bfs_observer3(int n) : preds(n), dists(n)
    {
    }
    
    void on_visit(int p, int v)
    {
        preds[v] = p;
        dists[v] = p >= 0 ? dists[p] + 1 : 0;
    }
    
    void on_discover(int p, int v) { }
    
    Array<int> preds;
    Array<int> dists;
};


void do_bfs(BFSIterator& bfs_it, int ns, const int *seeds, std::vector<int>& vs)
{    
    for (int i = 0; i < ns; ++i) 
        bfs_it.add_seed(seeds[i]);
            
    int v;
    while ((v = bfs_it.next()) >= 0)
    {
        vs.push_back(v+1);
    }        
}


void do_bfs(BFSIterator& bfs_it, int ns, const int *seeds, 
        std::vector<int>& vs, std::vector<int>& preds)
{
    for (int i = 0; i < ns; ++i)
        bfs_it.add_seed(seeds[i]);
    
    bfs_observer2 obs(bfs_it.num_nodes());
    
    int v;
    while ((v = bfs_it.next(obs)) >= 0)
    {
        vs.push_back(v+1);        
        preds.push_back(obs.preds[v]+1);
    }
}


void do_bfs(BFSIterator& bfs_it, int ns, const int *seeds, 
        std::vector<int>& vs, std::vector<int>& preds, std::vector<int>& dists)
{
    for (int i = 0; i < ns; ++i)
        bfs_it.add_seed(seeds[i]);
    
    bfs_observer3 obs(bfs_it.num_nodes());
    
    int v;
    while ((v = bfs_it.next(obs)) >= 0)
    {
        vs.push_back(v+1);        
        preds.push_back(obs.preds[v]+1);
        dists.push_back(obs.dists[v]);
    }
}



/***********************************************************
 *
 *  DFS 
 *
 **********************************************************/


struct dfs_observer2
{
    dfs_observer2(int n) : preds(n)
    {
    }
    
    void on_discover(int p, int v)
    {
        preds[v] = p;
    }
    
    void on_finish(int v) { }
    
    Array<int> preds;
};

struct dfs_observer3
{
    dfs_observer3(int n) : preds(n)
    {
        ford.reserve(n);
    }
    
    void on_discover(int p, int v)
    {
        preds[v] = p;
    }
    
    void on_finish(int v) 
    {
        ford.push_back(v);
    }
    
    Array<int> preds;
    std::vector<int> ford;
};

struct dfs_observer5
{
    dfs_observer5(int n) : time(0), preds(n), dtimes(n), ftimes(n)
    {
        ford.reserve(n);
    }
    
    void on_discover(int p, int v)
    {
        preds[v] = p;
        dtimes[v] = time++;
    }
    
    void on_finish(int v) 
    {
        ford.push_back(v);
        ftimes[v] = time++;
    }
    
    int time;
    
    Array<int> preds;
    std::vector<int> ford;
    
    Array<int> dtimes;
    Array<int> ftimes;
};




void do_dfs(DFSIterator& dfs_it, int ns, const int *seeds, std::vector<int>& vs)
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
                vs.push_back(v+1);
            }
        }
    }   
}

void do_dfs(DFSIterator& dfs_it, int ns, const int *seeds, std::vector<int>& vs,
        std::vector<int>& preds)
{        
    dfs_observer2 obs(dfs_it.num_nodes());
    
    for (int i = 0; i < ns; ++i) 
    {
        int s = seeds[i];
        if (!dfs_it.is_discovered(s))
        {
            dfs_it.set_seed(s);
            int v;
            while ((v = dfs_it.next(obs)) >= 0) 
            {
                vs.push_back(v+1);
                preds.push_back(obs.preds[v]+1);
            }
        }
    }   
}

void do_dfs(DFSIterator& dfs_it, int ns, const int *seeds, std::vector<int>& vs,
        std::vector<int>& preds, std::vector<int>& ford)
{        
    dfs_observer3 obs(dfs_it.num_nodes());
    
    for (int i = 0; i < ns; ++i) 
    {
        int s = seeds[i];
        if (!dfs_it.is_discovered(s))
        {
            dfs_it.set_seed(s);
            int v;
            while ((v = dfs_it.next(obs)) >= 0) 
            {
                vs.push_back(v+1);
                preds.push_back(obs.preds[v]+1);
            }
        }
    }  
    
    ford.swap(obs.ford);
    std::transform(ford.begin(), ford.end(), ford.begin(), 
            std::bind2nd(std::plus<int>(), 1));
}


void do_dfs(DFSIterator& dfs_it, int ns, const int *seeds, std::vector<int>& vs,
        std::vector<int>& preds, std::vector<int>& ford, 
        std::vector<int>& dtimes, std::vector<int>& ftimes)
{        
    dfs_observer5 obs(dfs_it.num_nodes());
    
    for (int i = 0; i < ns; ++i) 
    {
        int s = seeds[i];
        if (!dfs_it.is_discovered(s))
        {
            dfs_it.set_seed(s);
            int v;
            while ((v = dfs_it.next(obs)) >= 0) 
            {
                vs.push_back(v+1);
                preds.push_back(obs.preds[v]+1);                
            }
        }
    }  
    
    ford.swap(obs.ford);
    std::transform(ford.begin(), ford.end(), ford.begin(), 
            std::bind2nd(std::plus<int>(), 1));
    
    for (std::vector<int>::const_iterator it = vs.begin(); it != vs.end(); ++it)
    {
        int v = (*it) - 1;
        
        dtimes.push_back(obs.dtimes[v]);
        ftimes.push_back(obs.ftimes[v]);
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
        BFSIterator bfs_it(adjList);
        
        std::vector<int> preds;
        std::vector<int> dists;
        
        if (nlhs <= 1)
        {            
            do_bfs(bfs_it, ns, seeds, vs);        
        }    
        else if (nlhs == 2)
        {
            preds.reserve(G.nnodes());            
            
            do_bfs(bfs_it, ns, seeds, vs, preds); 
        }
        else if (nlhs == 3)
        {            
            preds.reserve(G.nnodes());            
            dists.reserve(G.nnodes());
            
            do_bfs(bfs_it, ns, seeds, vs, preds, dists);            
        }
        
        plhs[0] = v2mrow(vs);
        if (nlhs >= 2)
            plhs[1] = v2mrow(preds);
        if (nlhs >= 3)
            plhs[2] = v2mrow(dists);
        
    }
    else if (code == 'd')
    {
        DFSIterator dfs_it(adjList);
        
        std::vector<int> preds;
        std::vector<int> ford;
        std::vector<int> dtimes;
        std::vector<int> ftimes;     
                
        if (nlhs <= 1)
        {            
            do_dfs(dfs_it, ns, seeds, vs);
        }
        else if (nlhs == 2)
        {
            preds.reserve(G.nnodes());
            
            do_dfs(dfs_it, ns, seeds, vs, preds);
        }
        else if (nlhs == 3)
        {
            preds.reserve(G.nnodes());
            ford.reserve(G.nnodes());
            
            do_dfs(dfs_it, ns, seeds, vs, preds, ford);
        }
        else if (nlhs <= 5)
        {
            preds.reserve(G.nnodes());
            ford.reserve(G.nnodes());
            dtimes.reserve(G.nnodes());
            ftimes.reserve(G.nnodes());
            
            do_dfs(dfs_it, ns, seeds, vs, preds, ford, dtimes, ftimes);
        }
        
        plhs[0] = v2mrow(vs);
        if (nlhs >= 2)
            plhs[1] = v2mrow(preds);
        if (nlhs >= 3)
            plhs[2] = v2mrow(ford);
        if (nlhs >= 4)
            plhs[3] = v2mrow(dtimes);
        if (nlhs >= 5)
            plhs[4] = v2mrow(ftimes);
        
    }
}
   



