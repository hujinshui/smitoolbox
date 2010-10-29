/********************************************************************
 *
 *  array.h
 *
 *  The class to represent a fixed-length array
 *
 *  Created by Dahua Lin, on Oct 3, 2010
 *
 ********************************************************************/


#ifndef SMI_CLIB_ARRAY_H
#define SMI_CLIB_ARRAY_H

#include "commons.h"

namespace smi
{
 
template<typename T>
struct CRefMemory
{
    int n;
    const T *base;
    
    CRefMemory(size_t n_, const T *p) : n(n_), base(p) { }
};    
      
    
template<typename T>
struct RefMemory
{
    int n;
    T *base;
    
    RefMemory(size_t n_, T *p) : n(n_), base(p) { }
};
    

template<typename T>
RefMemory<T> crefmem(size_t n, const T *p)
{
    return RefMemory<T>(n, p);
}

template<typename T>
RefMemory<T> refmem(size_t n, T *p)
{
    return RefMemory<T>(n, p);
}
    

/************************************************
 *
 * Fixed-size Array classes
 *
 ************************************************/


template<typename T>
class CRefArray
{
public:
    SMI_DEFINE_COMMON_STL_TYPES(T)
    
public:
    CRefArray(size_type n, const value_type *data)
    : m_n(n), m_data(const_cast<T*>(data))
    {
    }
    
    size_type size() const
    {
        return m_n;
    }
    
    bool empty() const
    {
        return m_n == 0;
    }
    
    const_pointer data() const
    {
        return m_data;
    }
    
    const_reference operator[] (size_type i) const
    {
        return m_data[i];
    }
    
    const_pointer ptr(size_type i) const
    {
        return m_data + i;
    }
    
    const_iterator begin() const
    {
        return m_data;
    }
    
    const_iterator end() const
    {
        return m_data + m_n;
    }
    
    const_iterator rbegin() const
    {
        return m_data + (m_n - 1);
    }
    
    const_iterator rend() const
    {
        return m_data - 1;
    }
        
    const_reference front() const
    {
        return *m_data;
    }
    
    const_reference back() const
    {
        return m_data[m_n - 1];
    }        
    
public:
    
    bool elementwise_equal(const CRefArray<T>& rhs) const
    {
        if (m_n == rhs.m_n)
        {
            for (int i = 0; i < m_n; ++i) 
            {
                if (!(m_data[i] == rhs.m_data[i])) return false;
            }
            return true;
        }
        else
        {
            return false;
        }
    }
    
    bool bitwise_equal(const CRefArray<T>& rhs) const
    {
        if (m_n == rhs.m_n)
        {
            return ::memcmp(this->m_data, rhs.m_data, sizeof(value_type) * m_n) == 0;
        }
        else
        {
            return false;
        }        
    }
    
    bool reference_equal(const CRefArray<T>& rhs) const
    {
        return m_n == rhs.m_n && m_data == rhs.m_data;
    }    
    
protected:
    size_type m_n;
    value_type *m_data;
    
}; // end class CRefArray
    


template<typename T>
class RefArray : public CRefArray<T>
{
public:
    SMI_DEFINE_COMMON_STL_TYPES(T)
    
public:
    RefArray(size_type n, value_type *data) : CRefArray<T>(n, data)    
    {
    }
        
    const_pointer data() const
    {
        return this->m_data;
    }
    
    pointer data()
    {
        return this->m_data;
    }
    
    const_reference operator[] (size_type i) const
    {
        return this->m_data[i];
    }
    
    reference operator[] (size_type i) 
    {
        return this->m_data[i];
    }
    
    const_pointer ptr(size_type i) const
    {
        return this->m_data + i;
    }
    
    pointer ptr(size_type i) 
    {
        return this->m_data + i;
    }
    
    const_iterator begin() const
    {
        return this->m_data;
    }
    
    iterator begin() 
    {
        return this->m_data;
    }
    
    const_iterator end() const
    {
        return this->m_data + this->m_n;
    }
    
    iterator end()
    {
        return this->m_data + this->m_n;
    }
    
    const_iterator rbegin() const
    {
        return this->m_data + (this->m_n - 1);
    }
    
    iterator rbegin()
    {
        return this->m_data + (this->m_n - 1);
    }
    
    const_iterator rend() const
    {
        return this->m_data - 1;
    }
    
    iterator rend() 
    {
        return this->m_data - 1;
    }
        
    const_reference front() const
    {
        return *(this->m_data);
    }
    
    reference front() 
    {
        return *(this->m_data);
    }
    
    const_reference back() const
    {
        return this->m_data[this->m_n - 1];
    }        
    
    reference back()
    {
        return this->m_data[this->m_n - 1];
    }
    
public:
    
    void set_zeros()
    {
        ::memset(this->data(), 0, this->size() * sizeof(value_type));
    }
    
    void fill_value(const T& v)
    {
        pointer dst = this->data();
        size_type n = this->size();
        
        for (size_type i = 0; i < n; ++i) dst[i] = v;
    }
    
    template<typename TIter>
    void import_from(TIter it)
    {
        pointer dst = this->data();
        size_type n = this->size();
        
        for (size_type i = 0; i < n; ++i) dst[i] = *(it++);
    }
    
}; // end class RefArray



template<typename T>
class Array : public RefArray<T>
{
public:
    SMI_DEFINE_COMMON_STL_TYPES(T)
    
public:
    explicit Array(size_type n) : RefArray<T>(n, new T[n])
    {
    }
    
    explicit Array(size_type n, const T& v) : RefArray<T>(n, new T[n])
    {
        this->fill_value(v);
    }    
    
    ~Array()
    {
        delete[] this->m_data;
    }
    
    void swap(Array<T>& rhs)
    {
        size_type ts = this->m_n;
        this->m_n = rhs.m_n;
        rhs.m_n = ts;
        
        value_type *tp = this->m_data;
        this->m_data = rhs.m_data;
        rhs.m_data = tp;
    }
    
    
private:
    // disable copying
    
    Array(const Array<T>& );
    
    Array<T>& operator = (const Array<T>& );
    
}; // end class Array


typedef Array<bool> BoolArray;

}


#endif

