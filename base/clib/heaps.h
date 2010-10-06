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

#include "relation.h"
#include "data_structs.h"
#include "bintree.h"


namespace smi
{
    

/**
 * The class to represent a binary heap
 */    
template<typename TKey, typename HeapOrder=less<TKey> >    
class BinaryHeap
{
public:
    typedef TKey key_type;   
    typedef CompleteBinaryTree<int> btree_type;
    typedef typename btree_type::node_indicator trnode;
    
public:
    BinaryHeap(int cap) 
    : m_btree(cap), m_keys(cap), m_nodemap(cap), m_inheap(cap)    
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
    
    bool in_heap(int i) const
    {
        return i < m_keys.size() && m_inheap[i];
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
        return m_btree.root_value();
    }
    
    const btree_type& btree() const
    {
        return m_btree;
    }
    
public:
            
    void make_heap(int n, const key_type *src)
    {                        
        // copy key values and initialize the tree
        for (int i = 0; i < n; ++i) 
        {
            append_node(m_keys[i] = src[i]);            
        }
                                
        // heapify the tree
        if (n > 0)
        {
            int dir;
            trnode rend = m_btree.rev_end();
            
            for (trnode p = m_btree.last_nonleaf(); p != rend; --p)
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
        m_keys.clear();
    }
    
            
    void add_key(key_type v)
    {
        // append the key as last node
        append_node(v);
        
        // heapify the new node
        trnode p = m_btree.last_node();
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
            update_node(m_nodemap[i]);
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
                trnode pl = m_btree.last_node();
                swap_node(m_btree.root(), pl);
                
                // delete the last node
                m_inheap[m_btree[pl]] = false;
                m_btree.pop_node();
                
                // update the root node
                update_node(m_btree.root());
            }
            else
            {
                // delete the last node, which is also the root node
                m_inheap[m_btree.root_value()] = false; 
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
        trnode endp = m_btree.end();
        for (trnode p = m_btree.root(); p != endp; ++p)
        {
            dump_node_to_matlab(p);            
            mexPrintf(" ");
        }       
        mexPrintf("\n");
    }    

    void dump_node_to_matlab(trnode p) const
    {          
        double v = double(key_at_node(p));            
        mexPrintf("%g", v);
    }
    
    void dump_imap_to_matlab() const
    {
        for (int i = 0; i < m_keys.size(); ++i)
        {
            trnode p = m_nodemap[i];
            
            if (in_heap(i))
            {
                mexPrintf("(%d: %d -> %g) ", 
                        i, p.index(), double(key_at_node(p)));
            }
            else
            {
                mexPrintf("(%d: ~)", i); 
            }            
        }
        mexPrintf("\n");
    }

#endif
    
       
private:
    
    key_type key_at_node(trnode p) const
    {
        return m_keys[m_btree[p]];
    }
    
                      
    void append_node(key_type kv)
    {
        int new_index = m_keys.size();
        
        m_btree.push_node(new_index);
        m_nodemap[new_index] = m_btree.last_node();
        m_inheap[new_index] = true;
                
        m_keys.add(kv);        
    }
    
    
    void update_node(trnode p)
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
    
    bool superior_to_parent(trnode p) const
    {
        trnode par = p.parent();
        if (m_btree.has_node(par))
        {
            return _compare( key_at_node(p), key_at_node(par) );
        }
        else
        {
            return false;
        }                
    }
    
    int inferior_to_children(trnode p) const // return direction of the proper child
    {
        trnode lc = p.left_child();                     
        
        if (m_btree.has_node(lc))
        {
            key_type k = key_at_node(p);
            int dir = -1;
            
            key_type klc = key_at_node(lc);            
            if (_compare(klc, k))
            {
                k = klc;
                dir = 0;
            }
            
            trnode rc = lc.succ();
            if (m_btree.has_node(rc))
            {
                key_type krc = key_at_node(rc);
                if (_compare(krc, k))
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
    
    
    void upheap(trnode p) 
    {
        do 
        {
            p = swap_node(p, p.parent());
        } 
        while( superior_to_parent(p) );
    }
    
    void downheap(trnode p, int dir)
    {                        
        do
        {         
            p = swap_node(p, p.child(dir));    
        }
        while( (dir = inferior_to_children(p)) >= 0 );
    }
     
    trnode swap_node(trnode ps, trnode pt)
    {
        
#ifdef MATLAB_DUMP_HEAP_ACTION
        
        mexPrintf("swap ");
        dump_node_to_matlab(ps);
        mexPrintf(" <--> ");
        dump_node_to_matlab(pt);
        mexPrintf("\n");
        
#endif
        int si = m_btree[ps];
        int ti = m_btree[pt];
        
        m_btree[ps] = ti;
        m_btree[pt] = si;
        
        m_nodemap[si] = pt;
        m_nodemap[ti] = ps;        
               
        return pt;
    }           
            
private:
    HeapOrder _compare;  
    btree_type m_btree;  // binary tree of index: from node --> index
    
    SeqList<key_type> m_keys;   // map from index --> key
    Array<trnode> m_nodemap;   // map from index --> node 
    
    Array<bool> m_inheap;   // whether an index is in-heap
    
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

