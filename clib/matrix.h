/********************************************************************
 *
 *  vector.h
 *
 *  The header file for classes representing matrices
 *
 *  Created by Dahua Lin, on Sep 14, 2010
 *
 ********************************************************************/

#ifndef SMI_CLIB_MATRIX_H
#define SMI_CLIB_MATRIX_H

#include "vector.h"

namespace smi
{

template<typename T>
class MatrixCView
{
public:
    MatrixCView(const T *data, int m, int n)
    : m_data(const_cast<T*>(data)), m_nrows(m), m_ncols(n)
    {
    }
    
    int nrows() const
    {
        return m_nrows;
    }
    
    int ncols() const
    {
        return m_ncols;
    }
    
    int nelems() const
    {
        return m_nrows * m_ncols;
    }
    
    int first_nonsingleton_dim()
    {
        return m_nrows > 1 ? 0 : 1;
    }
    
    const T* data() const
    {
        return m_data;
    }
    
    const T& operator() (int i, int j) const
    {
        return m_data[i + j * m_nrows];
    }
    
    const T* ptr(int i, int j) const
    {
        return m_data + (i + j * m_nrows);
    }
    
    VectorCView<T> flat() const
    {
        return VectorCView<T>(m_data, nelems());
    }
        
    MultiVectorCView<T> columns() const
    {
        return MultiVectorCView<T>(m_data, m_ncols, m_nrows, m_nrows, 1);
    }
    
    MultiVectorCView<T> rows() const
    {
        return MultiVectorCView<T>(m_data, m_nrows, m_ncols, 1, m_nrows);
    }
    
protected:
    T *m_data;
    int m_nrows;
    int m_ncols;    
};
    

template<typename T>
class MatrixView : public MatrixCView<T>
{
public:    
    MatrixView(T *data, int m, int n)
    : MatrixCView<T>(data, m, n)
    {
    }
    
public:
    const T* data() const
    {
        return m_data;
    }
    
    T* data()
    {
        return m_data;
    }
    
    const T& operator() (int i, int j) const
    {
        return m_data[i + j * m_nrows];
    }
    
    T& operator() (int i, int j)
    {
        return m_data[i + j * m_nrows];
    }
    
    const T* ptr(int i, int j) const
    {
        return m_data + (i + j * m_nrows);
    }
    
    T* ptr(int i, int j)
    {
        return m_data + (i + j * m_nrows);
    }
    
    VectorCView<T> flat() const
    {
        return VectorCView<T>(m_data, nelems());
    }
    
    VectorView<T> flat() 
    {
        return VectorView<T>(m_data, nelems());
    }
        
    MultiVectorCView<T> columns() const
    {
        return MultiVectorCView<T>(m_data, m_ncols, m_nrows, m_nrows, 1);
    }
    
    MultiVectorView<T> columns()
    {
        return MultiVectorView<T>(m_data, m_ncols, m_nrows, m_nrows, 1);
    }
        
    MultiVectorCView<T> rows() const
    {
        return MultiVectorCView<T>(m_data, m_nrows, m_ncols, 1, m_nrows);
    }
    
    MultiVectorView<T> rows() 
    {
        return MultiVectorView<T>(m_data, m_nrows, m_ncols, 1, m_nrows);
    }
};
    
}

#endif
