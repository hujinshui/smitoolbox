/********************************************************************
 *
 *  marray.h
 *
 *  The class to represent a MATLAB array
 *
 *  Created by Dahua Lin, on Sep 14, 2010
 *
 ********************************************************************/

#ifndef SMI_CLIB_MARRAY_H
#define SMI_CLIB_MARRAY_H

#include <mex.h>

namespace smi
{
 
// The class to represent a MATLAB array    
class MArray 
{

public:
    MArray(const mxArray *A) : mxA(A)
    {        
    }
    
    int nrows() const
    {
        return (int)mxGetM(mxA);
    }
    
    int ncols() const
    {
        return (int)mxGetN(mxA);
    }
    
    int nelems() const
    {
        return (int)mxGetNumberOfElements(mxA);
    }
    
    int ndims() const
    {
        return (int)mxGetNumberOfDimensions(mxA);
    }
    
public:    
    bool is_empty() const
    {
        return mxIsEmpty(mxA);
    }
    
    bool is_matrix() const
    {
        return ndims() == 2;
    }
    
    bool is_vector() const
    {
        return ndims() == 2 && (nrows() == 1 || ncols() == 1);
    }
    
    bool is_sparse() const
    {
        return mxIsSparse(mxA);
    }
    
    bool is_complex() const
    {
        return mxIsComplex(mxA);
    }
    
    bool is_real_nonsparse_matrix() const
    {
        return is_matrix() && !is_complex() && !is_sparse();
    }
    
    bool is_scalar() const
    {
        return nelems() == 1 && !is_sparse();
    }
    
    bool is_real_double_scalar() const
    {
        return is_scalar() && is_double() && !is_complex();
    }
    
    bool is_logical_scalar() const 
    {
        return mxIsLogicalScalar(mxA);
    }
        
public:
    
    bool is_double() const
    {
        return mxIsDouble(mxA);
    }
    
    bool is_single() const
    {
        return mxIsSingle(mxA);
    }
    
    bool is_float() const
    {
        return is_double() || is_single();
    }
    
    bool is_numeric() const
    {
        return mxIsNumeric(mxA);
    }
    
    bool is_int8() const
    {
        return mxIsInt8(mxA);
    }
    
    bool is_int16() const
    {
        return mxIsInt16(mxA);
    }
    
    bool is_int32() const
    {
        return mxIsInt32(mxA);
    }
    
    bool is_int64() const
    {
        return mxIsInt64(mxA);
    }    
        
    bool is_uint8() const
    {
        return mxIsUint8(mxA);
    }
    
    bool is_uint16() const
    {
        return mxIsUint16(mxA);
    }
    
    bool is_uint32() const
    {
        return mxIsUint32(mxA);
    }
    
    bool is_uint64() const
    {
        return mxIsUint64(mxA);
    }
    
    bool is_char() const
    {
        return mxIsChar(mxA);
    }
    
    bool is_logical() const
    {
        return mxIsLogical(mxA);
    }
    
    bool is_struct() const
    {
        return mxIsStruct(mxA);
    }
    
    bool is_cell() const
    {
        return mxIsCell(mxA);
    }        
    
    
public:    
    double get_double_scalar() const
    {
        return mxGetScalar(mxA);
    }
    
    template<typename T>
    T get_scalar() const
    {
        return *(get_data<T>());
    }
    
    template<typename T>
    const T* get_data() const
    {
        return (const T*)mxGetData(mxA);
    }
    
    template<typename T>
    T* get_data() 
    {
        return (T*)mxGetData(mxA);
    }    
    
    template<typename T>
    VectorCView<T> to_vector() const
    {
        return VectorCView<T>(get_data<T>(), nelems());
    }
    
    template<typename T>
    VectorView<T> to_vector() 
    {
        return VectorView<T>(get_data<T>(), nelems());
    }
    
    template<typename T>
    MatrixCView<T> to_matrix() const
    {
        return MatrixCView<T>(get_data<T>(), nrows(), ncols());
    }
    
    template<typename T>
    MatrixView<T> to_matrix() 
    {
        return MatrixView<T>(get_data<T>(), nrows(), ncols());
    }
    
private:
    const mxArray *mxA;
    
}; // end class MArray
       
}

#endif

