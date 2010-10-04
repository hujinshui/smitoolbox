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

#include <string.h>

namespace smi
{
    
template<typename T>
class ConstArray
{
public:
    typedef T value_type;
    
public:
    ConstArray(int n, const value_type *data)
    : m_n(n), m_data(const_cast<T*>(data))
    {
    }
    
    int size() const
    {
        return m_n;
    }
    
    const value_type* data() const
    {
        return m_data;
    }
    
    const value_type& operator[] (int i) const
    {
        return m_data[i];
    }
    
    const value_type* ptr(int i) const
    {
        return m_data + i;
    }
        
protected:
    int m_n;
    value_type *m_data;
};
    

template<typename T>
class RefArray : public ConstArray<T>
{
public:
    typedef T value_type;
    
public:
    RefArray(int n, T *data) : ConstArray<T>(n, data)
    {
    }
    
    const value_type* data() const
    {
        return this->m_data;
    }
    
    value_type* data()
    {
        return this->m_data;
    }
    
    const value_type& operator[] (int i) const
    {
        return this->m_data[i];
    }
    
    value_type& operator[] (int i)
    {
        return this->m_data[i];
    }
    
    const value_type* ptr(int i) const
    {
        return this->m_data + i;
    }
    
    value_type* ptr(int i)
    {
        return this->m_data + i;
    }
    
    void set_zeros()
    {
        ::memset(this->m_data, 0, sizeof(value_type) * this->size());
    }      
    
    void fill_value(int v)
    {
        int n = this->size();
        value_type *a = this->data();
        for (int i = 0; i < n; ++i)
        {
            a[i] = v;
        }
    }
};


template<typename T>
class Array : public RefArray<T>
{
public:
    typedef T value_type;
    
public:
    explicit Array(int n) : RefArray<T>(n, new T[n])
    {
    }
    
    ~Array()
    {
        delete[] this->m_data;
    }
    
private:
    // disable copying
    
    Array(const Array<T>& );
    
    Array<T>& operator = (const Array<T>& );
    
};


template<typename T>
void memory_copy(const ConstArray<T>& src, RefArray<T>& dst)
{
    ::memcpy(dst.data(), src.data(), sizeof(T) * src.size());
}

}


#endif

