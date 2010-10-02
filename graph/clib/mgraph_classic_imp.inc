// The implementation of mgraph_classic

#include <vector>
#include <deque>

template<typename T>    
int smi::breadth_first_traverse(const smi::GNeighborHood<T>& G, int v0, int *r, bool *visited)
{
    std::deque<int> Q(G.nnodes());
    
    if (!visited[v0]) Q.push_back(v0);
    
    int l = 0;
    while (!Q.empty())
    {
        int s = Q.pop_front();        
        r[l++] = s;
        
        visited[s] = true;
        
        int n = G.neighbor_num(s);
        if (n > 0)
        {
            const int *nbs = G.neighbor_indices(s);
            
            for (int j = 0; j < n; ++j)
            {
                int t = nbs[j];
                if (!visited[t]) Q.push_back(t);                     
            }
        }                
    }
    
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
        
        dfs_entry(const smi::GNeighborHood<T>& G, int inode)
        {
            v = inode;
            n_nbs = G.neighbor_num();
            nbs = G.neighbor_indices();
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


template<typename T>
int smi::depth_first_traverse(const smi::GNeighborHood<T>& G, int v0, int *r, bool *visited)
{
    std::stack<smi_impl::dfs_entry> S(G.nnodes());
    
    int l = 0;    
    if (!visited[v0])
    {
        r[l++] = v0;
        S.push_back(dfs_entry(G, v0));
    }
    
    while (!S.empty())
    {
        smi_impl::dfs_entry& e = S.back();
        
        if (e.remain())
        {
            int v = e.next();
            if (!visited[v])
            {
                r[l++] = v;
                S.push_back(dfs_entry(G, v));
            }
        }
        else
        {
            S.pop_back();
        }
    }
}





