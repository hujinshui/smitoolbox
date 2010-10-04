/********************************************************************
 *
 *  heaps.h
 *
 *  The heap classes
 *
 *  Created by Dahua Lin, on Oct 4, 2010
 *
 ********************************************************************/


#ifndef SMI_CLIB_HEAPS_H
#define SMI_CLIB_HEAPS_H

#include "data_structs.h"


namespace smi
{

struct min_heap
{
    template<typename T>
    static bool compare(T x, T y) { return x < y; }
};


struct max_heap
{
    template<typename T>
    static bool compare(T x, T y) { return x > y; }
}
    
    

/**
 * The class to represent a binary heap
 */    
template<typename TKey, typename HeapOrder=min_heap>    
class BinaryHeap
{
public:
    typedef TKey key_type;
    typedef HeapOrder heap_order;
    
public:
    BinaryHeap(int cap, int n0, const TKey *init_keys)
    : m_capa(cap), m_n(n0)
    , m_keys(cap), m_stree(cap+1), m_imap(cap)
    {        
        for (int i = 0; i < n0; ++i) m_keys[i] = init_keys[i];
        
        make_heap();
    }
    
public:
    
    int capacity() const
    {
        return m_capa;
    }
    
    int size() const
    {
        return m_n;
    }
    
    key_type get_key(int i) const
    {
        return m_keys[i];
    }        
    
    int root_index() const
    {
        return m_stree[1];
    }
    
    key_type root_key() const
    {
        return get_key(root_index());
    }
    
    
public:
    
    int index_of_node(int p) const
    {
        return smap[p];
    }
        
    key_type key_of_node(int p) const
    {
        return get_key(index_of_node(p));
    }    
    
    int parent_node(int p) const
    {
        return p >> 1;
    }
    
    int left_child_node(int p) const
    {
        return p << 1;
    }
    
    int right_child_node(int p) const
    {
        return left_child_node(p) + 1;
    }
        
    int child_node(int p, int dir) const
    {
        int lc = left_child_node(p);
        return lc + (dir - 1);
    }
    
public:
    
    void add_key(key_type v)
    {
        m_keys[m_n] = v;        
        m_stree[m_n+1] = m_n;
        m_imap[m_n] = m_n+1;
        ++m_n;
        
        update_node(m_n);
    }
    
    
    void set_key(int i, key_type v) 
    {        
        if (m_keys[i] != v)
        {
            m_keys[i] = v;            
            update_node(m_imap[i]);
        }        
        
    }
    
    void delete_root()
    {
        if (m_n > 1)
        {
            swap_node(1, m_n);
            -- m_n;
            
            update_node(1);
        }
        else
        {
            m_n = 0;
        }
        
    }
    
    key_type extract_root()
    {
        key_type kv = root_key();
        delete_root();
        return kv;
    }
    
private:        
    
    void update_node(int p)
    {                
        if (superior_to_parent(p))
        {
            upheap(p);
        }
        else 
        {
            int dir = 0;
            if ((dir = inferior_to_children(p)) > 0)
            {
                downheap(p, dir);
            }
        }
    }
    
    bool superior_to_parent(int p) const
    {
        int par = parent_node(p);
        if (par)
        {
            return heap_order::compare( key_of_node(p), key_of_node(par) );
        }
        else
        {
            return false;
        }                
    }
    
    int inferior_to_children(int p) const
    {
        int lc = left_child_node(p);                       
        
        if (lc <= m_n)
        {
            k = key_of_node(p);
            dir = 0;
            
            int klc = key_of_node(lc);            
            if (heap_order::compare(klc, k))
            {
                k = klc;
                dir = 1;
            }
            
            rc = lc + 1;
            if (rc <= m_n)
            {
                int krc = key_of_node(rc);
                if (heap_order::compare(krc, k))
                {
                    k = krc;
                    dir = 2;
                }
            }
            
            return dir;                
        }
        else
        {
            return 0;
        }                
    }
    
    
    void upheap(int p) 
    {
        do 
        {
            p = swap_node(p, parent_of_node(p));
        } 
        while( superior_to_parent(p) );
    }
    
    void downheap(int p, int dir)
    {        
        do
        {
            p = swap_node(p, child_of_node(p, dir));
        }
        while( (dir = inferior_to_children(p)) > 0 );
    }
     
    int swap_node(int ps, int pt)
    {
        int si = index_of_node(ps);
        int ti = index_of_node(pt);
        
        m_stree[ps] = ti;
        m_stree[pt] = si;
        
        m_imap[si] = pt;
        m_imap[ti] = ps;
        
        return pt;
    }
    
    
    void make_heap()
    {        
        if (m_n > 0)
        {
            m_stree[0] = -1;
            for (int i = 0; i < m_n; ++i)
            {
                m_stree[i+1] = i;
                m_imap[i] = i+1;
            }
            
            int pl = parent_of_node(m_n);
            int dir = 0;
            for (int p = pl; p >= 1; --p)
            {
                if ( (dir = inferior_to_children(p)) > 0 )
                {
                    downheap(p, dir);
                }
            }
        }                
    }
    
    
private:
    BinaryHeap(const BinaryHeap<TKey>& );
    BinaryHeap<TKey>& operator = (const BinaryHeap<TKey>& );
    
    
private:
    int m_capa;
    int m_n;
    Array<key_type> m_keys;
    
    Array<int> m_stree; // use one-based index (for efficient computation of parent/children)
    Array<int> m_imap;
    
}; // end class BinaryHeap
    
    


}


#endif

