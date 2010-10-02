// The implementation of mgraph_classic

#include "mgraph_classic.h"
#include <string.h>

int smi::breadth_first_traverse(const smi::GNeighborHood& G, int v0, int *r, bool *visited)
{  
    int* Q = new int[G.nnodes()];
    int qb = 0;
    int qf = 0;

    if (!visited[v0])
    {
        Q[qb++] = v0;
        visited[v0] = true;
    }
    
    int l = 0;
    while (qf < qb)
    {
        int s = Q[qf++];
        r[l++] = s;
                
        int n = G.neighbor_num(s);
        if (n > 0)
        {
            const int *nbs = G.neighbor_nodes(s);
            
            for (int j = 0; j < n; ++j)
            {                
                int t = nbs[j];
                if (!visited[t])
                {
                    Q[qb++] = t;
                    visited[t] = true;
                }                    
            }
        }                
    }
    
    delete[] Q;
    
    return l;
}

namespace smi_impl
{
    struct dfs_entry
    {
        int v;          
        int n_nbs;
        const int *nbs;
        int c;
        
        dfs_entry()
        {
        }
        
        dfs_entry(const smi::GNeighborHood& G, int inode)
        {
            v = inode;
            n_nbs = G.neighbor_num(v);
            nbs = G.neighbor_nodes(v);
            c = 0;
        }
        
        bool remain() const
        {
            return c < n_nbs;
        }
        
        int next()
        {
            return nbs[c++];
        }
    };        
}


int smi::depth_first_traverse(const smi::GNeighborHood& G, int v0, int *r, bool *visited)
{
    smi_impl::dfs_entry *S = new smi_impl::dfs_entry[G.nnodes()];
    int sb = 0;
    
    int l = 0;    
    if (!visited[v0])
    {
        r[l++] = v0;
        visited[v0] = true;

        S[sb++] = smi_impl::dfs_entry(G, v0);
    }
    
    while (sb > 0)
    {
        smi_impl::dfs_entry& e = S[sb-1];

        if (e.remain())
        {
            int v = e.next();
            if (!visited[v])
            {
                r[l++] = v;
                visited[v] = true;
                
                S[sb++] = smi_impl::dfs_entry(G, v);
            }
        }
        else
        {
            -- sb;
        }
    }
    
    delete[] S;
    return l;
}


int smi::extract_connected_components(const smi::GNeighborHood& G, smi::GTraverseOrder gto, int *r, int *csiz)
{
    int n = G.nnodes();
    
    bool *visited = new int[n];
    ::memset(visited, 0, n * sizeof(bool));
    
    bool use_bf = (gto != smi::GTO_DF);
    
    int nc = 0;
    int *p = r;
    
    for (int i = 0; i < n; ++i)
    {
        if (!visited[i])
        {
            if (use_bf)
            {
                csiz[nc] = breadth_first_traverse(G, i, p, visited);
            }
            else
            {
                csiz[nc] = depth_first_traverse(G, i, p, visited);
            }      
            p += csiz[nc++];
        }                
    }
        
    delete[] visited;
    
    return nc;
}




