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
    
template<typename T> struct mttrait<T>;

template<> struct mttrait<int>
{
    const mxClassID class_id = mxINT32_CLASS;
    const bool is_float = false;
    const bool is_int = true;    
};

template<> struct mttrait<unsigned int>
{
    const mxClassID class_id = mxUINT32_CLASS;
    const bool is_float = false;
    const bool is_int = true;
};

template<> struct mttrait<float>
{
    const mxClassID class_id = mxSINGLE_CLASS;
    const bool is_float = true;
    const bool is_int = false;
};

template<> struct mttrait<double>
{
    const mxClassID class_id = mxDOUBLE_CLASS;
    const bool is_float = true;
    const bool is_int = false;
};

template<> struct mttrait<bool>
{
    const mxClassID class_id = mxLOGICAL_CLASS;
    const bool is_float = false;
    const bool is_int = false;
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




}


#endif


