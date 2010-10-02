// The implementation of mgraph_search

#include "mgraph_search.h"

int smi::gnode_trace_back(const int *prec, int v, int *path)
{    
    path[0] = v;
    int n = 1;
    
    int p;
    while ((p = prec[v]) >= 0)
    {
        path[n++] = v = p;        
    }       
    
    int l = 0;
    int r = n-1;
    while (l < r)
    {
        int tmp = path[l];
        path[l++] = path[r];
        path[r--] = tmp;        
    }
    
    return n;
}


/****************************
 *
 * BreadthFirstSearch
 *
 ****************************/

class smi::BreadthFirstSearch::Impl
{
public:
    Impl(GTraversal& gtr) : m_traversal(gtr)
    {
        m_queue = new int[gtr.neighborhood().nnodes()];
        m_qf = m_qb = 0;
    }
    
    ~Impl()
    {
        delete[] m_queue;
    }
    
    bool queue_empty() const
    {
        return m_qf == m_qb;
    }
    
    void add_start(int v)
    {
        if (!m_traversal.is_visited(v))
        {
            m_traversal.visit(-1, v);
        }
        m_queue[m_qb++] = v;
    }
    
    void enqueue(int prec, int v)
    {
        m_traversal.visit(prec, v);
        m_queue[m_qb++] = v;
    }
    
    int dequeue() 
    {
        return m_queue[m_qf++];
    }
    
private:
    GTraversal& m_traversal;
    int *m_queue;
    int m_qf;
    int m_qb;
};

void smi::BreadthFirstSearch::initialize(int ns, int *starts)
{
    _impl = new Impl(m_traversal);
    for (int i = 0; i < ns; ++i)
    {
        _impl->add_start(starts[i]);
    }
}

bool smi::BreadthFirstSearch::search(int vstop)
{
    if (vstop >= 0 && m_traversal.is_visited(vstop)) 
        return true;
    
    const GNeighborHood& G = m_traversal.NeighborHood();
    
    while (!_impl->empty_queue())
    {
        int s = _impl->dequeue();
        
        int n = G.neighbor_num(s);
        if (n > 0)
        {
            const int *nbs = G.neighbor_nodes(s);
            
            for (int j = 0; j < n; ++j)
            {                
                int t = nbs[j];
                if (!m_traversal.is_visited(t))
                {
                    _impl->enqueue(s, t);
                }                    
            }
        } 
        
        if (vstop >= 0 && m_traversal.is_visited(vstop)) 
            return true;
    }
    
    return false;
}


/****************************
 *
 * DepthFirstSearch
 *
 ****************************/

class smi::DepthFirstSearch::Impl
{
public:
    struct entry
    {
        int v;          
        int n_nbs;
        const int *nbs;
        int c;
        
        entry()
        {
        }
        
        entry(const smi::GNeighborHood& G, int inode)
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
    
public:
    Impl(GTraversal& gtr, int ns, const int *starts) 
    : m_traversal(gtr), m_G(gtr.neighborhood())
    {
        m_stack = new entry[gtr.neighborhood().nnodes()];
        m_sb = 0;
        
        m_ns = ns;
        m_starts = new int[m_ns];
        for (int i = 0; i < ns; ++i) m_starts[i] = starts[i];
                
        m_i0 = 0;
        next_stack();
    }
            
    ~Impl()
    {
        delete[] m_starts;
        delete[] m_stack;        
    }
    
    bool stack_empty()
    {
        return m_sb == 0;
    }
    
    void push(int prec, int v)
    {
        m_traversal.visit(prec, v);
        m_stack[m_sb++] = entry(m_G, v);
    }
    
    void pop()
    {
        -- m_sb;
        next_stack();
    }
    
    void next_stack()
    {
        while (m_i0 < m_ns)
        {
            int s = m_starts[m_i0++];
            if (!m_traversal.is_visited(s))
            {
                push(-1, s);
                break;
            }
        }
    }
    
    entry& top()
    {
        return m_stack[m_sb-1];
    }
        
private:
    GTraversal& m_traversal;
    GNeighorHood& m_G;
    
    int m_ns;
    int *m_starts;
    int m_i0;
    
    int *m_stack;
    int m_sb;
};

void smi::DepthFirstSearch::initialize(int ns, int *starts)
{
    _impl = new Impl(m_traversal, ns, starts);
}


void smi::DepthFirstSearch::search(int vstop)
{    
    if (vstop >= 0 && m_traversal.is_visited(vstop))            
        return true;
    
    while (sb > 0)
    {
        Impl::entry& e = _impl->top();

        if (e.remain())
        {
            int v = e.next();
            if (!m_traversal.visited[v])
            {
                _impl->push(e.v, v);
            }
        }
        else
        {
            _impl->pop();
        }
        
        if (vstop >= 0 && m_traversal.is_visited(vstop))
            return true;
    }
    
    return false;
}


