/********************************************************************
 *
 *  relation.h
 *
 *  Some structs/functors for relation judgment
 *
 *  Created by Dahua Lin, on Oct 5, 2010
 *
 ********************************************************************/

#ifndef SMI_CLIB_RELATION_H
#define SMI_CLIB_RELATION_H

namespace smi
{
    
    
template<typename T> 
struct less
{
    bool operator() (const T& x, const T& y) const
    {        
        return x < y;
    }
};


template<typename T> 
struct greater
{
    bool operator() (const T& x, const T& y) const
    {
        return x > y;
    }
};

  
template<typename T> 
struct less_equal
{
    bool operator() (const T& x, const T& y) const
    {
        return x <= y;
    }
};  
    

template<typename T> 
struct greater_equal
{
    bool operator() (const T& x, const T& y) const
    {
        return x >= y;
    }
}; 


template<typename T>
struct equal
{
    bool operator() (const T& x, const T& y) const
    {
        return x == y;
    }
};


template<typename T>
struct not_equal
{
    bool operator() (const T& x, const T& y) const
    {
        return x != y;
    }
};

    
    
};


#endif

