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

#include <vector>

namespace smi
{    
        
    
/**
 * The class to represent a complete binary tree
 */     
template<typename T>    
class CompleteBinaryTree    
{    
public:
    
    typedef std::vector<T> container_type;
        
    typedef T value_type;
<<<<<<< local
    typedef container_type::size_type size_type;
    typedef container_type::const_reference const_reference;
    typedef container_type::reference reference;
    typedef container_type::const_iterator const_iterator;
    typedef container_type::iterator iterator;
        
=======
    typedef T& reference;
    typedef const T& const_reference;
    typedef size_t size_type;
    
>>>>>>> other
public:
    
    class node_indicator
    {        
    public:                
        node_indicator() { }
        
        node_indicator(size_type p) : _p(p) { }
        
        size_type index() const 
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
            return (_p << 1) + (size_type)dir;
        }
        
        node_indicator& operator++()
        {
            ++_p;
            return *this;
        }
        
        node_indicator operator++(int)
        {
            return _p++;            
        }
        
        node_indicator& operator--() 
        {
            --_p;
            return *this;
        }
        
        node_indicator operator--(int)
        {
            return _p--;
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
        size_type _p;        
<<<<<<< local
=======
    };
>>>>>>> other
        
        
    }; // end class node_indicator
                
    
public:
<<<<<<< local
    CompleteBinaryTree() : m_stree(1)
    {        
    }    
    
    CompleteBinaryTree(size_type initcap)
=======
    CompleteBinaryTree(size_type cap)
    : m_capa(cap), m_n(0), m_stree(cap+1)
>>>>>>> other
    {
        m_stree.reserve(initcap + 1);
        m_stree.push(value_type());
    }    
    
    size_type capacity() const
    {
        return m_stree.capacity() - 1;
    }
    
<<<<<<< local
    size_type size() const
=======
    size_type capacity() const
>>>>>>> other
    {
<<<<<<< local
        return m_stree.size() - 1;
=======
        return m_capa;
    }
    
    size_type size() const
    {
        return m_n;
>>>>>>> other
    }
    
    bool empty() const
    {
        return m_stree.size() <= 1;
    }       
            
    node_indicator root_node() const
    {
        return 1;
    }    
           
    node_indicator last_node() const
    {
        return size();
    }        
    
    node_indicator last_nonleaf() const
    {
        return last_node().parent();
    }
        
    
public:
    
    
    bool has_node(node_indicator p) const
    {
        return p.index() > 0 && p.index() <= m_n;
    }
    
    const_reference operator[] (node_indicator p) const
    {
        return m_stree[p.index()];
    }
    
    reference operator[] (node_indicator p) 
    {
        return m_stree[p.index()];
    }
    
    const_iterator get_iter(node_indicator p) const
    {
        return m_stree.begin() + p.index();
    }
    
    iterator get_iter(node_indicator p)
    {
        return m_stree.begin() + p.index();
    }
            
    
    const_iterator begin() const
    {
        return m_stree.begin() + 1;
    }
    
    iterator begin()
    {
        return m_stree.begin() + 1;
    }
    
    const_iterator end() const
    {
        return m_stree.end();
    }
    
    iterator end()
    {
        return m_stree.end();
    }            
    
    const_reverse_iterator rbegin() const
    {
        return m_stree.rbegin();
    }
    
    reverse_iterator rbegin() 
    {
        return m_stree.rbegin();
    }
        
    const_reverse_iterator rend() const
    {
        return m_stree.rend() - 1;
    }
    
    const_reverse_iterator rend()
    {
        return m_stree.rend() - 1;
    }
        
public:
        
    /**
     * add a node as the last node in the tree
     */
    void push_node(const value_type& e)
    {        
        m_stree.push_back(e);
    }
    
    
    /**
     * remove the last node
     */
    void pop_node()
    {
        m_stree.pop_back();
    }
             
    
    /**
     * clear all nodes
     */
    void clear()
    {
        m_stree.clear();
        m_stree.push_back(value_type());
    }
    
    
private:
<<<<<<< local
    container_type m_stree;  
=======
    size_type m_capa;
    size_type m_n;
>>>>>>> other
    
}; // end class CompleteBinaryTree


}

#endif

