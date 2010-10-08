/********************************************************************
 *
 *  disjointset.h
 *
 *  The class to implement disjoint-set structure
 *
 *  Created by Dahua Lin, on Oct 8, 2010
 *
 ********************************************************************/


#ifndef SMI_CLIB_DISJOINT_SET_H
#define SMI_CLIB_DISJOINT_SET_H


#include <vector>

namespace smi
{

    
class DisjointSetForest
{
public:
        
    typedef std::vector<int>::size_type size_type;
    typedef size_type index_type;
    
    struct entry
    {
        entry() { }
        entry(index_type p) : parent(p), rank(0) { }
        
        index_type parent;
        size_type rank;
    };
    
public:
    DisjointSetForest() 
    { 
    }    
    
    explicit DisjointSetForest(size_type initcapa)    
    {
        m_entries.reserve(initcapa);
    }           
    
public:
    
    index_type parent(index_type x) const
    {
        return m_entries[x].parent;
    }
        
    bool is_root(index_type x) const
    {
        return parent(x) == x;
    }            
    
    index_type find_root(index_type x) const
    {
        index_type p = parent(x);
        while (p != x)
        {
            p = parent(p);
        }
        return p;
    }    
        
    bool in_same_set(size_type x, size_type y) const
    {
        return find_root(x) == find_root(y);
    }
    
    
    void add_singletons(size_type n)
    {
        index_type idx = m_entries.size();
        for (size_type i = 0; i < n; ++i)
        {
            m_entries.push_back(entry(idx++));
        }
    } 
    
    bool merge(size_type x, size_type y)
    {
        index_type xroot = find_root(x);
        index_type yroot = find_root(y);
        
        if (xroot != yroot)
        {
            size_type xrank = rank(xroot);
            size_type yrank = rank(yroot);
            
            if (xrank < yrank)
            {
                m_entries[xroot].parent = yroot;
            }
            else if (xrank > yrank)
            {
                m_entries[yroot].parent = xroot;
            }
            else
            {
                m_entries[yroot].parent = xroot;
                m_entries[xroot].rank = xrank + 1;
            }
        }
        else
        {
            return false;
        }
    }
    
private:
    size_type rank(index_type x)
    {
        return m_entries[x].rank;
    }
    
    
private:
    std::vector<entry> m_entries;
    
}; // end class DisjointSetForest
    
    

}

#endif

