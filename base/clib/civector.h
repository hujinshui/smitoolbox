/********************************************************************
 *
 *  civector.h
 *
 *  The consistently indexing container
 *
 *  Created by Dahua Lin, on Oct 7, 2010
 *
 ********************************************************************/


#ifndef SMI_CLIB_CIVECTOR_H
#define SMI_CLIB_CIVECTOR_H

#include <vector>

namespace smi
{

/**
 * A dynamic array with consistent indexing
 *
 * Each element has identified with a unique integer index, 
 * which will not change even when the internal data-structure
 * is re-organized.
 *
 */    
template<typename T>    
class CIVector    
{
public:
    
    struct entry_t
    {
        T value;
        bool in;
        
        entry_t() { }
        
        entry_t(T v) : value(v), in(true) { }
    };
        
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    
    typedef typename std::vector<entry_t>::size_type size_type;
    typedef size_type index_type;
    
public:
    
    CIVector() : m_size(0) { }
    
    CIVector(size_type initcapa) : m_size(0) 
    {
        m_entries.reserve(initcapa);
    }
    
    void reserve(size_type capa)
    {
        m_entries.reserve(capa);
    }    
    
    size_type capacity() const
    {
        return m_entries.capacity();
    }
        
    size_type size() const 
    {
        return m_size;
    }
                
    bool empty() const
    {
        return m_size == 0;
    }
    
    bool has_index(index_type i) const
    {
        return i < m_entries.size() && m_entries[i].in;
    }
    
    index_type next_index() const
    {
        return m_entries.size();
    }
        
    const_reference operator[] (index_type i) const
    {
        return m_entries[i].value;
    }
        
    reference operator[] (index_type i)
    {
        return m_entries[i].value;
    }    
    
    
    index_type add(const T& v)
    {
        m_entries.push_back(entry_t(v));
        ++ m_size;
    }
    
    void remove(index_type i)  // pre-condition: this->has_index(i)
    {
        m_entries[i].in = false;
        -- m_size;
    }
    
    void clear()
    {
        m_entries.clear();
        m_size = 0;
    }    
        
private:    
    std::vector<entry_t> m_entries; 
    size_type m_size;
    
}; // end class CIVector
    

};

#endif

