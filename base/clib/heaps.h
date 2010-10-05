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


#ifdef MATLAB_INSPECT_HEAPS
#include <mex.h>
#endif

// to dump heap actions, further define the macro MATLAB_DUMP_HEAP_ACTION

#include "bintree.h"

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
};
    
    

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
    BinaryHeap(int cap) : m_btree(cap), m_keys(cap)    
    {        
    }
    
public:
    
    int capacity() const
    {
        return m_btree.capacity();
    }
    
    int size() const
    {
        return m_btree.size();
    }
    
    bool empty() const
    {
        return m_btree.empty();
    }
        
    key_type get_key(int i) const
    {
        return m_keys[i];
    }        
        
    key_type root_key() const
    {
        return get_key(root_index());
    }            
    
    int root_index() const
    {
        return m_btree.root_index();
    }
    
    const AssoBinaryTree& btree() const
    {
        return m_btree;
    }
    
public:
            
    void make_heap(int n, const key_type *src)
    {                        
        // copy key values
        for (int i = 0; i < n; ++i) m_keys[i] = src[i];
        
        // reset the tree
        m_btree.make_tree(n);
                
        // heapify the tree
        if (n > 0)
        {
            int dir;
            for (int p = m_btree.last_nonleaf(); p > 0; --p)
            {
                if ((dir = inferior_to_children(p)) >= 0)
                {
                    downheap(p, dir);
                }
            }            
        }
    }
    
    
    void clear()
    {
        m_btree.clear();        
    }
    
            
    void add_key(key_type v)
    {
        // add to key array
        m_keys[size()] = v;
        
        // push a new node to the tree
        m_btree.push_node();
        
        // heapify the new node
        int p = m_btree.last_node();
        if (superior_to_parent(p))
        {
            upheap(p);
        }
    }
    
    
    void set_key(int i, key_type v) 
    {        
        if (m_keys[i] != v)
        {
            // update the key value
            m_keys[i] = v;       
            
            // heapify the updated node
            update_node(m_btree.node_of_index(i));
        }                
    }
    
    void delete_root()
    {
        if (!m_btree.empty())
        {
            int n = m_btree.size();
            if (n > 1)
            {
                // swap root with last node, and update
                swap_node(m_btree.root(), m_btree.last_node());
                
                // delete the last node
                m_btree.pop_node();
                
                // update the root node
                update_node(m_btree.root());
            }
            else
            {
                // delete the last node, which is also the root node
                m_btree.pop_node();
            }                        
        }
    }
    
    key_type extract_root(int& i)
    {
        i = root_index();
        key_type kv = get_key(i);
        
        delete_root();
        
        return kv;
    }    
            
        
#ifdef MATLAB_INSPECT_HEAPS

    void dump_heap_to_matlab() const
    {   
        for (int p = m_btree.root(); p <= m_btree.last_node(); ++p)
        {
            dump_node_to_matlab(p);            
            mexPrintf(" ");
        }       
        mexPrintf("\n");
    }    

    void dump_node_to_matlab(int p) const
    {          
        double v = double(key_at_node(p));            
        mexPrintf("%g", v);
    }
    
    void dump_imap_to_matlab() const
    {
        for (int i = 0; i < size(); ++i)
        {
            int p = m_btree.node_of_index(i);
            mexPrintf("(%d -> %g) ", p, double(key_at_node(p)));
        }
        mexPrintf("\n");
    }

#endif
    
       
private:        
    
    key_type key_at_node(int p) const
    {
        return get_key(m_btree.index_at_node(p));
    }    
    
    void update_node(int p)
    {                
        int dir;
        if (superior_to_parent(p))
        {
            upheap(p);
        }
        else if ((dir = inferior_to_children(p)) >= 0)
        {
            downheap(p, dir);
        }
    }
    
    bool superior_to_parent(int p) const
    {
        int par = m_btree.parent(p);
        if (par)
        {
            return heap_order::compare( key_at_node(p), key_at_node(par) );
        }
        else
        {
            return false;
        }                
    }
    
    int inferior_to_children(int p) const // return direction of the proper child
    {
        int maxp = m_btree.last_node();
        int lc = m_btree.left_child(p);                     
        
        if (lc <= maxp)
        {
            key_type k = key_at_node(p);
            int dir = -1;
            
            key_type klc = key_at_node(lc);            
            if (heap_order::compare(klc, k))
            {
                k = klc;
                dir = 0;
            }
            
            int rc = lc + 1;
            if (rc <= maxp)
            {
                key_type krc = key_at_node(rc);
                if (heap_order::compare(krc, k))
                {
                    k = krc;
                    dir = 1;
                }
            }
            
            return dir;                
        }
        else
        {
            return -1;
        }                
    }
    
    
    void upheap(int p) 
    {
        do 
        {
            p = swap_node(p, m_btree.parent(p));
        } 
        while( superior_to_parent(p) );
    }
    
    void downheap(int p, int dir)
    {                        
        do
        {         
            p = swap_node(p, m_btree.child(p, dir));    
        }
        while( (dir = inferior_to_children(p)) >= 0 );
    }
     
    int swap_node(int ps, int pt)
    {
        
#ifdef MATLAB_DUMP_HEAP_ACTION
        
        mexPrintf("swap ");
        dump_node_to_matlab(ps);
        mexPrintf(" <--> ");
        dump_node_to_matlab(pt);
        mexPrintf("\n");
        
#endif
        
        m_btree.swap_node(ps, pt);        
        return pt;
    }
            
private:
    AssoBinaryTree m_btree;
    Array<key_type> m_keys;        
    
}; // end class BinaryHeap
    
    
/**
 * Extract the root-values sequentially from a well-formed heap
 *
 * @param heap the heap from which the keys are to be extracted
 * @param n the number of values to be extracted (n <= heap.size())
 * @param keys the buffer to store the output keys (can be NULL)
 * @param indices the buffer to store the output indices (can be NULL)
 *
 * @remarks A heap-sort can be accomplished by make_heap -> extract_all
 */
template<typename THeap>
void heap_extract_n(THeap& heap, int n, typename THeap::key_type *keys, int *indices = 0)
{
    int i = 0;
    
    if (keys != 0)
    {
        if (indices != 0)        
            for (int t = 0; t < n; ++t) 
                keys[t] = heap.extract_root(indices[t]);
        else
            for (int t = 0; t < n; ++t) 
                keys[t] = heap.extract_root(i);
    }
    else
    {
        if (indices != 0)        
            for (int t = 0; t < n; ++t) 
                heap.extract_root(indices[t]);
        else
            for (int t = 0; t < n; ++t) 
                heap.extract_root(i);
    }
}



}


#endif

