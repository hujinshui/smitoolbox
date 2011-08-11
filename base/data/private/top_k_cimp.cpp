/********************************************************************
 *
 *  top_k_cimp.cpp
 *
 *  The C++ mex implementation of top_k.m
 *
 *  Created by Dahua Lin, on Nov 11, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

#include <functional>
#include <algorithm>


using namespace bcs;
using namespace bcs::matlab;


template<typename T>
struct indexed_value
{
    T v;
    int i;
    
    void set(int idx, T val)
    {
        v = val;
        i = idx;
    }
    
    bool operator < (const indexed_value<T>& rhs) const
    {
        return v < rhs.v;
    }
    
    bool operator > (const indexed_value<T>& rhs) const
    {
        return v > rhs.v;
    }
    
    bool operator == (const indexed_value<T>& rhs) const
    {
        return v == rhs.v;
    }
};


template<typename T>
inline void export_with_inds(const T *src, indexed_value<T>* dst, index_t n)
{    
    for (index_t i = 0; i < n; ++i)
    {
        dst[i].set(i, src[i]);
    }
}


template<typename T, class Comp>
inline void top_k(aview1d<T> temp, aview1d<T> r, Comp comp)
{
    T *ws = temp.pbase();
    index_t n = temp.nelems();
    index_t K = r.nelems();
            
    std::nth_element(ws, ws + (K-1), ws + n, comp);
    std::sort(ws, ws + K, comp);
    
    for (index_t i = 0; i < K; ++i)
    {
        r(i) = temp(i);
    }
}


template<typename T, class Comp>
inline void top_k(aview1d<indexed_value<T> > temp, 
        aview1d<T> R, aview1d<double> I, Comp comp)
{
    indexed_value<T> *ws = temp.pbase();
    index_t n = temp.nelems();
    index_t K = R.nelems();
            
    std::nth_element(ws, ws + (K-1), ws + n, comp);
    std::sort(ws, ws + K, comp);
    
    for (index_t i = 0; i < K; ++i)
    {
        R(i) = temp(i).v;
        I(i) = temp(i).i + 1;
    }
}



template<typename T, template<typename U> class TComp>
inline marray do_top_k(const_marray mX, size_t K)
{    
    index_t m = (index_t)mX.nrows();
    index_t n = (index_t)mX.ncolumns();
        
    if (n == 1)
    {
        array1d<T> temp(m);
        
        marray mR = create_marray<T>(K,1);
        aview1d<T> R = view1d<T>(mR);
        
        copy(view1d<T>(mX), temp);
        top_k((aview1d<T>)temp, R, TComp<T>());
        
        return mR;
    }
    else
    {
        marray mR = create_marray<T>(K, n);
        array1d<T> temp(m);
        
        caview2d<T, column_major_t> X = view2d<T>(mX);
        aview2d<T, column_major_t> R = view2d<T>(mR);
        
        for (index_t i = 0; i < n; ++i)
        {
            copy(X.column(i), temp);
            top_k((aview1d<T>)temp, R.column(i), TComp<T>());
        }
        
        return mR;
    }
}


template<typename T, template<typename U> class TComp>
inline void do_indexed_top_k(const_marray mX, index_t K, marray& mR, marray& mI)
{    
    index_t m = (index_t)mX.nrows();
    index_t n = (index_t)mX.ncolumns();
    
    if (n == 1)
    {                
        mR = create_marray<T>((size_t)K, 1);
        mI = create_marray<double>((size_t)K, 1);
        array1d<indexed_value<T> > temp(m);

        caview1d<T> X = view1d<T>(mX);
        aview1d<T> R = view1d<T>(mR);
        aview1d<double> I = view1d<double>(mI);

        export_with_inds(X.pbase(), temp.pbase(), m);
        top_k((aview1d<indexed_value<T> >)temp, 
                R, I, TComp<indexed_value<T> >() );
    }
    else
    {
        mR = create_marray<T>((size_t)K, (size_t)n);  
        mI = create_marray<double>((size_t)K, (size_t)n);
        array1d<indexed_value<T> > temp(m);
        
        caview2d<T, column_major_t> X = view2d<T>(mX);
        aview2d<T, column_major_t> R = view2d<T>(mR);
        aview2d<double, column_major_t> I = view2d<double>(mI);
        
        for (index_t i = 0; i < n; ++i)
        {
            export_with_inds(&(X(0, i)), temp.pbase(), m);
            top_k((aview1d<indexed_value<T> >)temp, 
                    R.column(i), I.column(i), TComp<indexed_value<T> >() );
        } 
    }    
}



template<template<typename U> class TComp>
inline marray do_top_k_dtype(const_marray mX, index_t K)
{    
    switch (mX.class_id())
    {
        case mxDOUBLE_CLASS:
            return do_top_k<double, TComp>(mX, K);
        
        case mxSINGLE_CLASS:
            return do_top_k<float, TComp>(mX, K);
            
        case mxINT32_CLASS:
            return do_top_k<int32_t, TComp>(mX, K);
            
        case mxUINT32_CLASS:
            return do_top_k<uint32_t, TComp>(mX, K);
            
        default:
            throw mexception("top_k:invalidarg",
                "The X should be of type double, single, int32, or uint32.");
    }    
}

template<template<typename U> class TComp>
inline void do_indexed_top_k_dtype(const_marray mX, index_t K, marray& mR, marray& mI)
{    
    switch (mX.class_id())
    {
        case mxDOUBLE_CLASS:
            do_indexed_top_k<double, TComp>(mX, K, mR, mI);
            break;
        
        case mxSINGLE_CLASS:
            do_indexed_top_k<float, TComp>(mX, K, mR, mI);
            break;
            
        case mxINT32_CLASS:
            do_indexed_top_k<int32_t, TComp>(mX, K, mR, mI);
            break;
            
        case mxUINT32_CLASS:
            do_indexed_top_k<uint32_t, TComp>(mX, K, mR, mI);
            break;
        
        default:
            throw mexception("top_k:invalidarg",
                "The X should be of type double, single, int32, or uint32.");
    }    
}



/**
 * main entry:
 *
 * [0]:  X:     the input data matrix
 * [1]:  K:     the number of top elements
 * [2]:  code:  1 - max, 2 - min
 * 
 * output
 * [0]:  R:     top values
 * [1]:  I:     top indices (optional)
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const_marray mX(prhs[0]);
    const_marray mK(prhs[1]);
    const_marray mCode(prhs[2]);    
    
    index_t K = (index_t)mK.get_scalar<double>();
    int code = (int)mCode.get_scalar<double>();
    
    marray mR, mI;
    
    if (nlhs <= 1)
    {
        if (code == 1)
        {
            mR = do_top_k_dtype<std::greater>(mX, K);
        }
        else
        {
            mR = do_top_k_dtype<std::less>(mX, K);
        }
        
        plhs[0] = mR.mx_ptr();
    }
    else
    {
        if (code == 1)
        {
            do_indexed_top_k_dtype<std::greater>(mX, K, mR, mI);
        }
        else 
        {
            do_indexed_top_k_dtype<std::less>(mX, K, mR, mI);
        }
        
        
        plhs[0] = mR.mx_ptr();
        plhs[1] = mI.mx_ptr();
    }
               
}


BCSMEX_MAINDEF
      
        
