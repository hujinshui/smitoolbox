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
mxArray* create_matlab_matrix(int m, int n)
{
    return mxCreateNumericMatrix(m, n, mttrait<T>::class_id, mxREAL);
}

template<>
mxArray* create_matlab_matrix<bool>(int m, int n)
{
    return mxCreateLogicalMatrix(m, n);
}

template<typename T>
mxArray* create_matlab_scalar(T v)
{
    mxArray* M = create_matlab_matrix<T>(1, 1);
    *((T*)mxGetData(M)) = v;
    return M;
}


}


#endif


