/********************************************************************
 *
 *  heaps.h
 *
 *  The heap classes
 *
 *  Created by Dahua Lin, on Oct 4, 2010
 *
 *  Note: to inspect the behavior of a binary heap, one can 
 *  pre-define the macro SMI_BINARY_HEAP_INSPECTION
 *
 ********************************************************************/


#ifndef SMI_CLIB_HEAPS_H
#define SMI_CLIB_HEAPS_H

#include "civector.h"
#include "bintree.h"
#include <functional>


namespace smi
{
    

/**
 * The class to represent a binary heap
 */    
template<typename TKey, typename TOrd=std::less<TKey> >
class BinaryHeap
{    
       
public:    
    typedef TKey key_type;   
    typedef TOrd key_order;
        
    typedef CIVector<int>::index_type index_type;
    
    typedef CompleteBinaryTree<index_type> btree_type;
    typedef typename btree_type::size_type size_type;
    typedef typename btree_type::node_indicator trnode;
    
    struct entry
    {
        key_type key;
        trnode node;
        
        entry() { }        
        entry(const key_type& k, trnode nd) : key(k), node(nd) { }
    };
        
    typedef CIVector<entry> entry_container_type;
    
    
#ifdef SMI_BINARY_HEAP_INSPECTION
    
public:
    class IMonitor
    {
    public:
        virtual ~IMonitor() { }
        
        virtual void on_node_swaped(const BinaryHeap& H, index_type i, index_type j) = 0;
        
        virtual void on_node_appended(const BinaryHeap& H, index_type i) = 0;
        
        virtual void on_upheaping(const BinaryHeap& H, index_type i) = 0;
        
        virtual void on_downheaping(const BinaryHeap& H, index_type i) = 0;                
        
        virtual void on_key_setting(const BinaryHeap& H, index_type i, const key_type& kv) = 0;
        
        virtual void on_key_adding(const BinaryHeap& H, const key_type& kv) = 0;
        
        virtual void on_heap_making(const BinaryHeap& H) = 0;
        
        virtual void on_heap_made(const BinaryHeap& H) = 0;
        
        virtual void on_clearing(const BinaryHeap& H) = 0;
        
        virtual void on_cleared(const BinaryHeap& H) = 0;
        
        virtual void on_root_deleting(const BinaryHeap& H) = 0;
        
        virtual void on_root_deleted(const BinaryHeap& H) = 0;
        
    }; // end class Monitor
    
    
    void set_monitor(IMonitor *mon)
    {
        m_mon = mon;
    }
    
private:
    IMonitor* m_mon;

#endif
    
        
    
public:
    BinaryHeap()
    {
        #ifdef SMI_BINARY_HEAP_INSPECTION
                m_mon = 0;
        #endif
    }    
    
    BinaryHeap(size_type cap) 
    {
        m_btree.reserve(cap);
        m_map.reserve(cap);
        
        #ifdef SMI_BINARY_HEAP_INSPECTION
                m_mon = 0;
        #endif
    }
    
public:
    
    size_type capacity() const
    {
        return m_btree.capacity();
    }
    
    size_type size() const
    {
        return m_btree.size();
    }
    
    bool empty() const
    {
        return m_btree.empty();
    }
        
    const key_type& get_key(index_type i) const
    {
        return m_map[i].key;
    }        
    
    bool is_inheap(index_type i) const
    {
        return m_map.has_index(i);
    }
            
    const key_type& root_key() const
    {
        return key_at_node(m_btree.root_node());
    }            
     
    
    const index_type root_index() const
    {
        return m_btree[m_btree.root_node()];
    }
    
    const btree_type& btree() const
    {
        return m_btree;
    }
    
    const key_type& key_at_node(trnode p) const
    {        
        return get_key(m_btree[p]);
    }
    
    trnode get_node_by_index(index_type i) const
    {
        return m_map[i].node;
    }
    
public:
            
    template<typename InputIter>
    void make_heap(InputIter it_begin, InputIter it_end)
    {                        
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_heap_making(*this);
        #endif
        
        // copy key values and initialize the tree
        for (InputIter it = it_begin; it != it_end; ++it) 
        {
            append_node(*it);
        }
                                
        // heapify the tree
        
        trnode rend = m_btree.pre_root_node();            
        for (trnode p = m_btree.last_nonleaf() ; p != rend; --p)
        {
            int dir;
            if ((dir = inferior_to_children(p)) >= 0)
            {
                downheap(p, dir);
            }
        }   
        
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_heap_made(*this);
        #endif
    }
    
    
    void clear()
    {
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_heap_clearing(*this);
        #endif
        
        m_btree.clear(); 
        m_map.clear();
        
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_heap_cleared(*this);
        #endif
    }
    
            
    index_type add_key(const key_type& v)
    {
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_key_adding(*this, v);
        #endif
        
        // append the key as last node
        index_type idx = append_node(v);
        
        // heapify the new node
        trnode p = m_btree.last_node();
        if (superior_to_parent(p))
        {
            upheap(p);
        }
        
        return idx;
    }
    
    
    void set_key(int i, const key_type& v) 
    {        
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_key_setting(*this, i, v);
        #endif
        
        if (m_map[i].key != v)
        {                        
            // update the key value
            m_map[i].key = v;      
                                    
            // heapify the updated node
            update_node(m_map[i].node);            
        }                
    }
    
    void delete_root()
    {
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_root_deleting(*this);
        #endif
        
        if (!m_btree.empty())
        {                        
            index_type ridx = root_index();
            
            int n = m_btree.size();            
            if (n > 1)
            {                                
                // swap root with last node, and update
                trnode pl = m_btree.last_node();
                swap_node(m_btree.root_node(), pl);
                
                // delete the last node
                m_btree.pop_node();
                
                // update the root node
                update_node(m_btree.root_node());
            }
            else
            {
                // delete the last node, which is also the root node
                m_btree.pop_node();
            } 
            
            m_map.remove(ridx);
            
            #ifdef SMI_BINARY_HEAP_INSPECTION
                    if (m_mon != 0)
                        m_mon->on_root_deleted(*this);
            #endif
        }
    }
                           
       
private:                
                      
    index_type append_node(key_type kv)
    {
        index_type new_index = m_map.next_index();
        
        m_btree.push_node(new_index);
        m_map.add(entry(kv, m_btree.last_node()));
                    
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_node_appended(*this, new_index);
        #endif
                
        return new_index;
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
        if (!par.is_null())
        {
            return m_ord( key_at_node(p), key_at_node(par) );
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
            if (m_ord(klc, k))
            {
                k = klc;
                dir = 0;
            }
            
            trnode rc = lc.succ();
            if (m_btree.has_node(rc))
            {
                key_type krc = key_at_node(rc);
                if (m_ord(krc, k))
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
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_upheaping(*this, m_btree[p]);
        #endif
        
        do 
        {
            p = swap_node(p, p.parent());
        } 
        while( superior_to_parent(p) );
    }
    
    void downheap(trnode p, int dir)
    {      
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_downheaping(*this, m_btree[p]);
        #endif
                
        do
        {         
            p = swap_node(p, p.child(dir));    
        }
        while( (dir = inferior_to_children(p)) >= 0 );
    }
     
    trnode swap_node(trnode ps, trnode pt)
    {        
        index_type si = m_btree[ps];
        index_type ti = m_btree[pt];
        
        m_btree[ps] = ti;
        m_btree[pt] = si;
        
        m_map[si].node = pt;
        m_map[ti].node = ps;  
        
        #ifdef SMI_BINARY_HEAP_INSPECTION
                if (m_mon != 0)
                    m_mon->on_node_swaped(*this, si, ti);
        #endif
               
        return pt;
    }           
            
private:
    key_order m_ord;  
    btree_type m_btree;         // node -> index    
    entry_container_type m_map;   // index -> (key, node)       
    
}; // end class BinaryHeap
    
    

}


#endif

