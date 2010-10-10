/********************************************************************
 *
 *  functors.h
 *
 *  The header file provides a series of useful functors
 *
 *  Created by Dahua Lin, on Oct 10, 2010
 *
 ********************************************************************/


#ifndef SMI_CLIB_FUNCTORS_H
#define SMI_CLIB_FUNCTORS_H

#include <functional>

// Convenient functors


template<typename TArray>
struct array_based_mapper
{
    typedef int argument_type;
    typedef typename TArray::value_type result_type;
    
    const TArray& _array;
    
    array_based_mapper(const TArray& arr) : _array(arr) { }
    
    result_type operator() (argument_type i) const { return _array[i]; }
};


template<typename TArray>
array_based_mapper<TArray> arr_map(const TArray& arr)
{
    return array_based_mapper<TArray>(arr);
}


template<typename TFunctor1, typename TFunctor2>
struct unary_chain_functor2
{
    typedef TFunctor1 first_func_type;
    typedef TFunctor2 second_func_type;
    
    typedef typename first_func_type::argument_type argument_type;
    typedef typename second_func_type::result_type result_type;
            
    first_func_type first_func;
    second_func_type second_func;    
    
    unary_chain_functor2(const first_func_type& f1, const second_func_type& f2)
    : first_func(f1), second_func(f2)
    {
    }
     
    result_type operator() (argument_type x)
    {
        return f2(f1(x));
    }
};

template<typename TFunctor1, typename TFunctor2, typename TFunctor3>
struct unary_chain_functor3
{
    typedef TFunctor1 first_func_type;
    typedef TFunctor2 second_func_type;
    typedef TFunctor3 third_func_type;
    
    typedef typename first_func_type::argument_type argument_type;
    typedef typename third_func_type::result_type result_type;
            
    first_func_type first_func;
    second_func_type second_func;    
    third_func_type third_func;
    
    unary_chain_functor3(
            const first_func_type& f1, 
            const second_func_type& f2,
            const third_func_type& f3)
    : first_func(f1), second_func(f2), third_func(f3)
    {
    }
     
    result_type operator() (argument_type x)
    {
        return f3(f2(f1(x)));
    }
};



template<typename TFunctor1, typename TFunctor2>
inline unary_chain_functor2<TFunctor1, TFunctor2> unary_chain(
        const TFunctor1 &f1, 
        const TFunctor2 &f2)
{
    return unary_chain_functor2<TFunctor1, TFunctor2>(f1, f2);
}


template<typename TFunctor1, typename TFunctor2, typename TFunctor3>
inline unary_chain_functor3<TFunctor1, TFunctor2, TFunctor3> unary_chain(
        const TFunctor1 &f1, 
        const TFunctor2 &f2,
        const TFunctor3 &f3)
{
    return unary_chain_functor3<TFunctor1, TFunctor2, TFunctor3>(f1, f2, f3);
}


#endif



        