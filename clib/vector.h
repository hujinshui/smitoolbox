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

namespace smi
{
    
template<typename T>
class VectorCView
{
public:
    VectorCView(const T* data, int n) : m_data(const_cast<T*>(data)), m_n(n), m_intv(1)
    {
    }
    
    VectorCView(const T* data, int n, int intv) : m_data(const_cast<T*>(data)), m_n(n), m_intv(intv)
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
    
    const T& operator() (int i) const
    {
        return m_data[i * m_intv];
    }
    
    const T* ptr(int i) const
    {
        return m_data + i * m_intv;
    }
       
protected:
    T *m_data;  // the pointer to data
    int m_n;    // the number of elements
    int m_intv; // the interval between elements    
};
    

template<typename T>
class VectorView : public VectorCView<T>
{
public:
    VectorView(T *data, int n) : VectorCView<T>(data, n)
    {
    }
    
    VectorView(T *data, int n, int intv) : VectorCView<T>(data, n, intv)
    {
    }
    
    const T* data() const
    {
        return this->m_data;
    }
    
    T* data()
    {
        return this->m_data;
    }
    
    const T& operator() (int i) const
    {
        return this->m_data[i * this->m_intv];
    }
    
    T& operator() (int i) 
    {
        return this->m_data[i * this->m_intv];
    }
    
    const T* ptr(int i) const
    {
        return this->m_data + i * this->m_intv;
    }
    
    T* ptr(int i)
    {
        return this->m_data + i * this->m_intv;
    }
};


template<typename T>
class MultiVectorCView
{
public:
    MultiVectorCView(const T *data, int n, int len)
    : m_data(const_cast<T*>(data)), m_num(n), m_vlen(len), m_vintv(len), m_eintv(1)
    {
    }
    
    MultiVectorCView(const T *data, int n, int len, int vintv, int eintv)
    : m_data(const_cast<T*>(data)), m_num(n), m_vlen(len), m_vintv(vintv), m_eintv(eintv)
    {
    }
    
    int nvecs() const
    {
        return m_num;
    }
    
    int veclen() const
    {
        return m_vlen;
    }
    
    int vec_interval() const
    {
        return m_vintv;
    }
    
    int elem_interval() const
    {
        return m_eintv;
    }
    
    VectorCView<T> operator[] (int i) const
    {
        return VectorCView<T>(m_data + i * m_vintv, m_vlen, m_eintv);
    }
            
protected:
    T *m_data;
    int m_num;
    int m_vlen;
    int m_vintv;
    int m_eintv;    
};


template<typename T>
class MultiVectorView : public MultiVectorCView<T>
{
public:
    MultiVectorView(T *data, int n, int len)
    : MultiVectorCView<T>(data, n, len)
    {
    }
    
    MultiVectorView(T *data, int n, int len, int vintv, int eintv)
    : MultiVectorCView<T>(data, n, len, vintv, eintv)
    {
    }
    
    VectorCView<T> operator[] (int i) const
    {
        return VectorCView<T>(this->m_data + i * this->m_vintv, this->m_vlen, this->m_eintv);
    }
    
    VectorView<T> operator[] (int i) 
    {
        return VectorView<T>(this->m_data + i * this->m_vintv, this->m_vlen, this->m_eintv);
    }
};


template<typename T>
MultiVectorCView<T> emulate_multi(const VectorCView<T>& v, int n)
{
    return MultiVectorCView<T>(v.data(), n, v.nelems(), 0, v.interval());
}

template<typename T>
MultiVectorView<T> emulate_multi(VectorView<T>& v, int n) 
{
    return MultiVectorView<T>(v.data(), n, v.nelems(), 0, v.interval());
}

}


#endif

