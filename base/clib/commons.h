/********************************************************************
 *
 *  commons.h
 *
 *  The header file offers common facilities for SMI C++ library
 *
 *  Created by Dahua Lin, on Oct 6, 2010
 *
 ********************************************************************/

#ifndef SMI_CLI8_COMMONS_H
#define SMI_CLIB_COMMONS_H

#include <stddef.h>
#include <string.h>

#include <stdexcept>

#define SMI_DEFINE_COMMON_STL_TYPES(T) \
    typedef T value_type; \
    typedef T& reference;  \
    typedef const T& const_reference;  \
    typedef T* pointer;  \
    typedef const T* const_pointer;  \
    typedef size_t size_type;  \
    typedef ptrdiff_t difference_type;  \
    typedef pointer iterator;  \
    typedef const_pointer const_iterator; 


namespace smi
{

typedef unsigned char byte;   

};

#endif

