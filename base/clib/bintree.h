/********************************************************************
 *
 *  bintree.h
 *
 *  The class to represent a binary tree
 *
 *  Created by Dahua Lin, on Oct 4, 2010
 *
 ********************************************************************/

#ifndef SMI_CLIB_BINTREE_H
#define SMI_CLIB_BINTREE_H

#include "array.h"

namespace smi
{    
       
/**
 * The class to represent an associative binary tree
 *
 * The tree structure is maintained in an array.
 * Each node is attached an index, which can be used to 
 * refer to associative information in external world.
 *
 * The nodes: 1, ..., n
 * The indice: 0, ..., n-1
 */     
class AssoBinaryTree    
{
public:
    AssoBinaryTree(int cap)
    : m_capa(cap), m_n(0), m_stree(cap+1), m_imap(cap)
    {
    }
    
    int capacity() const
    {
        return m_capa;
    }
    
    int size() const
    {
        return m_n;
    }
    
    bool empty() const
    {
        return m_n == 0;
    }
    
public:
    
    int root_index() const
    {
        return m_stree[root()];
    }
    
    int index_at_node(int p) const
    {
        return m_stree[p];
    }
    
    int node_of_index(int i) const
    {
        return m_imap[i];
    }
    
    int root() const
    {
        return 1;
    }
    
    int parent(int p) const
    {
        return p >> 1; 
    }
    
    int left_child(int p) const
    {
        return p << 1;
    }
    
    int right_child(int p) const
    {
        return (p << 1) + 1;
    }
    
    int child(int p, int dir) const // dir: 0 - left, 1 - right
    {
        return (p << 1) + dir;
    }
    
    int last_node() const
    {
        return m_n;
    }        
    
    int last_nonleaf() const
    {
        return parent(last_node());
    }
    
public:
    
    /**
     * make a tree (arrange nodes at its natural order)
     */
    void make_tree(int n)
    {
        m_n = n;
        m_stree[0] = -1;    // this node is not used
        
        for (int i = 0; i < n; ++i)
        {
            m_stree[i+1] = i;
            m_imap[i] = i+1;
        }
    }
    
    /**
     * add a node as the last node in the tree
     */
    void push_node()
    {
        m_stree[m_n+1] = m_n;
        m_imap[m_n] = m_n + 1;
        ++m_n;
    }
    
    
    /**
     * remove the last node
     */
    void pop_node()
    {
        -- m_n;
    }
    
    /**
     * swap the position of two nodes
     */
    void swap_node(int p1, int p2)
    {
        int i1 = index_at_node(p1);
        int i2 = index_at_node(p2);
        
        m_stree[p1] = i2;
        m_stree[p2] = i1;
        
        m_imap[i1] = p2;
        m_imap[i2] = p1;
    }        
    
    
    void clear()
    {
        m_n = 0;
    }
    
private:
    int m_capa;
    int m_n;
    
    Array<int> m_stree; // map: node --> index
    Array<int> m_imap;  // map: index --> node
};


}

#endif

