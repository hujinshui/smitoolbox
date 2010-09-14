/********************************************************************
 *
 *  vector.h
 *
 *  The header file for classes representing vectors
 *
 *  Created by Dahua Lin, on Sep 14, 2010
 *
 ********************************************************************/

#ifndef SMI_CLIB_VECTORS_H
#define SMI_CLIB_VECTORS_H

#include <mex.h>

namespace smi
{
    
template<typename T>
class ConstVector
{
public:
    ConstVector(const T* data, int n) : m_data(const_cast<T*>(data)), m_n(n), m_intv(1)
    {
    }
    
    ConstVector(const T* data, int n, int intv) : m_data(const_cast<T*>(data)), m_n(n), m_intv(intv)
    {
    }
        
    int nelems() const
    {
        return m_n;
    }
    
    int interval() const
    {
        return m_intv;
    }
    
    const T* data() const
    {
        return m_data;
    }
    
    const T& operator[] (int i) const
    {
        return m_data[i * m_intv];
    }
    
    const T* ptr(int i) const
    {
        return m_data + i * m_intv;
    }
    
    mxArray* to_matlab_row() const
    {
        
    }
    
    mxArray* to_matlab_column() const
    {
    }
    
    
protected:
    T *m_data;  // the pointer to data
    int m_n;    // the number of elements
    int m_intv; // the interval between elements    
};
    

template<typename T>
class Vector : public ConstVector<T>
{
public:
    Vector(int n) : ConstVector<T>(new T[n], n), m_own(true)
    {
    }
    
    Vector(T *data, int n) : ConstVector<T>(data, n), m_own(false)
    {
    }
    
    Vector(T *data, int n, int intv) : ConstVector<T>(data, n, intv), m_own(false)
    {
    }

    ~Vector()
    {
        if (m_own)
        {
            delete[] m_data;
        }        
    }
    
    const T* data() const
    {
        return m_data;
    }
    
    T* data()
    {
        return m_data;
    }
    
    const T& operator[] (int i) const
    {
        return m_data[i * m_intv];
    }
    
    T& operator[] (int i) 
    {
        return m_data[i * m_intv];
    }
    
    const T* ptr(int i) const
    {
        return m_data + i * m_intv;
    }
    
    T* ptr(int i)
    {
        return m_data + i * m_intv;
    }
    
private:
    bool m_own;     // whether it owns the memory
}






}

#endif

