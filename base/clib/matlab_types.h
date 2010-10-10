/********************************************************************
 *
 *  matlab_types.h
 *
 *  The facilities for dealing with MATLAB types
 *
 *  Created by Dahua Lin, on Sep 14, 2010
 *
 ********************************************************************/

#ifndef SMI_CLIB_MATLAB_TYPES
#define SMI_CLIB_MATLAB_TYPES

#include <mex.h>
#include <string.h>
#include <iterator>

namespace smi
{
    
template<typename T> struct mttrait;

template<> struct mttrait<int>
{
    static const mxClassID class_id = mxINT32_CLASS;
    static const bool is_float = false;
    static const bool is_int = true;    
};

template<> struct mttrait<unsigned int>
{
    static const mxClassID class_id = mxUINT32_CLASS;
    static const bool is_float = false;
    static const bool is_int = true;
};

template<> struct mttrait<float>
{
    static const mxClassID class_id = mxSINGLE_CLASS;
    static const bool is_float = true;
    static const bool is_int = false;
};

template<> struct mttrait<double>
{
    static const mxClassID class_id = mxDOUBLE_CLASS;
    static const bool is_float = true;
    static const bool is_int = false;
};

template<> struct mttrait<bool>
{
    static const mxClassID class_id = mxLOGICAL_CLASS;
    static const bool is_float = false;
    static const bool is_int = false;
};


template<typename T>
inline mxArray* create_matlab_matrix(int m, int n)
{
    return mxCreateNumericMatrix(m, n, mttrait<T>::class_id, mxREAL);
}

inline mxArray* create_empty_matrix()
{
    return mxCreateDoubleMatrix(0, 0, mxREAL);
}

template<>
inline mxArray* create_matlab_matrix<bool>(int m, int n)
{
    return mxCreateLogicalMatrix(m, n);
}


template<typename T>
inline mxArray* src_to_matlab_matrix(int m, int n, const T *src)
{
    mxArray *mx = create_matlab_matrix<T>(m, n);
    T *dst = (T*)mxGetData(mx);    
    ::memcpy(dst, src, m * n * sizeof(T));
    return mx;
}

template<typename T>
inline mxArray* src_to_matlab_matrix_cond(int m, int n, const T *src)
{
    if (src != 0)
    {
        return src_to_matlab_matrix<T>(m, n, src);
    }
    else
    {
        return create_matlab_matrix<T>(0, 0);
    }
}



template<typename T>
inline mxArray* create_matlab_scalar(T v)
{
    mxArray* M = create_matlab_matrix<T>(1, 1);
    *((T*)mxGetData(M)) = v;
    return M;
}


template<typename TIter>
inline mxArray* iter_to_matlab_row(TIter it, int n)
{
    typedef typename std::iterator_traits<TIter>::value_type T;
    
    mxArray *mx = create_matlab_matrix<T>(1, n);
    T *dst = (T*)mxGetData(mx);
    for (int i = 0; i < n; ++i) dst[i] = *(it++);    
    return mx;
}

template<typename TIter, typename TGen>
inline mxArray* iter_to_matlab_row(TIter it, int n, TGen g)
{
    typedef typename TGen::result_type T;
        
    mxArray *mx = create_matlab_matrix<T>(1, n);
    T *dst = (T*)mxGetData(mx);
    for (int i = 0; i < n; ++i) dst[i] = g(*(it++));    
    return mx;
}



template<typename T>
struct offset_to_index
{
    typedef int argument_type;
    typedef T result_type;
    
    result_type operator() (argument_type i) const { return i+1; }
};



}


#endif


