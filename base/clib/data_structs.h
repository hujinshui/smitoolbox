/********************************************************************
 *
 *  data_struct.h
 *
 *  Template classes of some basic data structure.
 *
 *  Created by Dahua Lin, on Oct 03, 2010
 *
 ********************************************************************/


#ifndef SMI_CLIB_DATA_STRUCTS_H
#define SMI_CLIB_DATA_STRUCTS_H

#include <string.h>
#include <stdexcept>

namespace smi
{
    



class BoolArray
{
public:
    typedef int size_type;
    
public:
    explicit BoolArray(size_type n) : m_n(n), m_vals(new bool[n])
    {
    }
    
    BoolArray(size_type n, bool v) : m_n(n), m_vals(new bool[n])
    {
        fill_value(v);
    }
    
    ~BoolArray()
    {
        delete[] m_vals;
    }
    
    size_type size() const
    {
        return m_n;
    }
    
    int nbytes() const
    {
        return size() * sizeof(bool);
    }
    
    bool empty() const
    {
        return m_n == 0;
    }
    
    bool operator[] (int i) const
    {
        return m_vals[i];
    }
    
    bool& operator[] (int i) 
    {
        return m_vals[i];
    }
    
public:
    
    void fill_value(bool v)
    {
        if (!v) 
        {
            ::memset(m_vals, 0, nbytes());
        }
        else
        {
            int n = size();
            for (int i = 0; i < n; ++i) m_vals[i] = true;
        }
    }
    
    void copy_from(const bool *src)
    {
        ::memcpy(m_vals, src, nbytes());
    }
    
    bool operator == (const BoolArray& rhs) const
    {
        return m_n == rhs.m_n 
                && ::memcmp(m_vals, rhs.m_vals, nbytes()) == 0;
    }
    
    bool operator != (const BoolArray& rhs) const
    {
        return !(operator == (rhs));
    }
    
    bool all_true() const
    {
        for (int i = 0; i < m_n; ++i) 
        {
            if (!m_vals[i]) return false;
        }
        return true;
    }
    
    bool any_true() const
    {
        for (int i = 0; i < m_n; ++i)
        {
            if (m_vals[i]) return true;
        }
        return false;
    }
    
    bool all_false() const
    {
        return !any_true();
    }
    
    bool any_false() const
    {
        return !all_true();
    }
    
private:
    BoolArray(const BoolArray& );
    BoolArray& operator = (const BoolArray& );
    
private:
    bool *m_vals;
    size_type m_n;
    
}; // end class BoolArray


/**
 * The class serve as a fixed-capacity container base for stack or queue
 */
template<typename T>
class SQBase
{
public:
    typedef T value_type;
    typedef size_t size_type;
    
    typedef const T& const_reference;
    typedef T& reference;
    
public:
    SQBase(int c) : m_data(new T[c]), m_ifront(0), m_iback(-1), m_capa(c)
    {
    }
    
    ~SQBase()
    {
        delete[] m_data;
    }
    
    size_type capacity() const
    {
        return m_capa;
    }
    
    bool empty() const
    {
        return m_iback < m_ifront;
    }
    
    size_type size() const
    {
        int s = m_iback - m_ifront + 1;
        return s > 0 ? size_type(s) : 0;
    }
    
    reference front()
    {
        return m_data[m_ifront];
    }
    
    const_reference front() const
    {
        return m_data[m_ifront];
    }
    
    reference back()
    {
        return m_data[m_iback];
    }
    
    const_reference back() const
    {
        return m_data[m_iback];
    }
    
    void pop_front()
    {
        ++ m_ifront;
    }
    
    void pop_back()
    {
        -- m_iback;
    }
    
    void push_back(const T& e)
    {        
        m_data[++m_iback] = e;
    }                       
    
    void clear()
    {
        m_ifront = 0;
        m_iback = -1;
    }
    
private:
    SQBase(const SQBase<T>& );
    SQBase<T>& operator = (const SQBase<T>& );
    
private:
    T *m_data;
    int m_ifront;
    int m_iback;
    size_type m_capa;
    
}; // end class SQBase


    
}


#endif



