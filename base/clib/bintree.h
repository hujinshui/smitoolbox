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
 * The class to represent a complete binary tree
 */     
template<typename T>    
class CompleteBinaryTree    
{
public:
    typedef T value_type;
    
public:
    
    class node_indicator
    {
    public:
        node_indicator() { }
        
        node_indicator(int p) : _p(p) { }
        
        int index() const 
        { 
            return _p; 
        }        
        
        node_indicator pred() const
        {
            return _p - 1;
        }
        
        node_indicator succ() const
        {
            return _p + 1;
        }
        
        node_indicator parent() const
        {
            return _p >> 1;
        }
        
        node_indicator left_child() const
        {
            return _p << 1;
        }
        
        node_indicator right_child() const
        {
            return (_p << 1) + 1;
        }
        
        node_indicator child(int dir) const
        {
            return (_p << 1) + dir;
        }
        
        node_indicator& operator++() // prefix
        {
            ++_p;
            return *this;
        }
        
        node_indicator& operator--() // prefix
        {
            --_p;
            return *this;
        }
        
        bool operator == (node_indicator rhs) const
        {
            return _p == rhs._p;
        }
        
        bool operator != (node_indicator rhs) const
        {
            return _p != rhs._p;
        }
                
    private:
        int _p;        
    };
        
    
    
public:
    CompleteBinaryTree(int cap)
    : m_capa(cap), m_n(0), m_stree(cap+1)
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
    
    const value_type& operator[] (node_indicator p) const
    {
        return m_stree[p.index()];
    }
    
    value_type& operator[] (node_indicator p) 
    {
        return m_stree[p.index()];
    }
    
    bool has_node(node_indicator p) const
    {
        return p.index() > 0 && p.index() <= m_n;
    }    
    
    node_indicator root() const
    {
        return 1;
    }    
    
    const value_type& root_value() const
    {
        return m_stree[1];
    }
    
    
    node_indicator last_node() const
    {
        return m_n;
    }        
    
    node_indicator last_nonleaf() const
    {
        return last_node().parent();
    }
    
    
    node_indicator rev_end() const
    {
        return 0;
    }    
    
    node_indicator end() const
    {
        return m_n + 1;
    }
    
    
public:
        
    /**
     * add a node as the last node in the tree
     */
    void push_node(const value_type& e)
    {        
        m_stree[++m_n] = e;
    }
    
    
    /**
     * remove the last node
     */
    void pop_node()
    {
        -- m_n;
    }
             
    
    /**
     * clear all nodes
     */
    void clear()
    {
        m_n = 0;
    }
    
    
private:
    int m_capa;
    int m_n;
    
    Array<value_type> m_stree;    
};


}

#endif

