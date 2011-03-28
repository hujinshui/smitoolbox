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


template<typename T, class TIndexer>
inline void export_with_inds(
        const const_aview1d<T, TIndexer>& src, 
        aview1d<indexed_value<T> >& dst)
{
    int n = (int)src.nelems();
    
    for (int i = 0; i < n; ++i)
    {
        dst(i).set(i, src(i));
    }
}


template<typename T, class TIndexer, class Comp>
inline void top_k(aview1d<T> temp, aview1d<T, TIndexer> r, Comp comp)
{
    T *ws = temp.pbase();
    size_t n = temp.nelems();
    size_t K = r.nelems();
            
    std::nth_element(ws, ws + (K-1), ws + n, comp);
    std::sort(ws, ws + K, comp);
    
    for (index_t i = 0; i < (index_t)K; ++i)
    {
        r(i) = temp(i);
    }
}


template<typename T, class TIndexer, class Comp>
inline void top_k(aview1d<indexed_value<T> > temp, 
        aview1d<T, TIndexer> R, aview1d<double, TIndexer> I, Comp comp)
{
    indexed_value<T> *ws = temp.pbase();
    size_t n = temp.nelems();
    size_t K = R.nelems();
            
    std::nth_element(ws, ws + (K-1), ws + n, comp);
    std::sort(ws, ws + K, comp);
    
    for (index_t i = 0; i < (index_t)K; ++i)
    {
        R(i) = temp(i).v;
        I(i) = temp(i).i + 1;
    }
}



template<typename T, template<typename U> class TComp>
inline marray vec_top_k(const_marray mX, size_t K, int adim)
{    
    size_t n = mX.nelems();
    array1d<T> temp(n);
    
    marray mR = (adim == 1 ? create_marray<T>(K,1) : create_marray<T>(1,K));    
    aview1d<T> R = view1d<T>(mR);
    
    temp << view1d<T>(mX);
    top_k(temp, R, TComp<T>());
    
    return mR;
}

template<typename T, template<typename U> class TComp>
inline std::pair<marray, marray> vec_indexed_top_k(const_marray mX, size_t K, int adim)
{    
    size_t n = mX.nelems();
    array1d<indexed_value<T> > temp(n);
    
    marray mR, mI;
    
    if (adim == 1)
    {
        mR = create_marray<T>(K, 1);
        mI = create_marray<double>(K, 1);
    }
    else
    {
        mR = create_marray<T>(1, K);
        mI = create_marray<double>(1, K);
    }
    
    aview1d<T> R = view1d<T>(mR);
    aview1d<double> I = view1d<double>(mI);
    
    export_with_inds(view1d<T>(mX), temp);
    top_k(temp, R, I, TComp<indexed_value<T> >() );
    
    return std::make_pair(mR, mI);
}



template<typename T, template<typename U> class TComp>
marray mat_top_k(const_marray mX, size_t K, int adim)
{        
    size_t m = mX.nrows();
    size_t n = mX.ncolumns();
    
    if (adim == 1)  // columns
    {
        marray mR = create_marray<T>(K, n);         
        array1d<T> temp(m);
        
        const_aview2d<T, column_major_t> X = view2d<T>(mX);
        aview2d<T, column_major_t> R = view2d<T>(mR);
        
        for (size_t i = 0; i < n; ++i)
        {
            temp << X.column(i);
            top_k(temp, R.column(i), TComp<T>());
        }
        
        return mR;        
    }
    else  // rows
    {
        marray mR = create_marray<T>(m, K);
        array1d<T> temp(n);
        
        const_aview2d<T, column_major_t> X = view2d<T>(mX);
        aview2d<T, column_major_t> R = view2d<T>(mR);
        
        for (size_t i = 0; i < m; ++i)
        {
            temp << X.row(i);
            top_k(temp, R.row(i), TComp<T>());
        }
        
        return mR;
    }             
}


template<typename T, template<typename U> class TComp>
std::pair<marray, marray> mat_indexed_top_k(const_marray mX, size_t K, int adim)
{    
    size_t m = mX.nrows();
    size_t n = mX.ncolumns();
    
    if (adim == 1)  // columns
    {
        marray mR = create_marray<T>(K, n);  
        marray mI = create_marray<double>(K, n);
        array1d<indexed_value<T> > temp(m);
        
        const_aview2d<T, column_major_t> X = view2d<T>(mX);
        aview2d<T, column_major_t> R = view2d<T>(mR);
        aview2d<double, column_major_t> I = view2d<double>(mI);
        
        for (size_t i = 0; i < n; ++i)
        {
            export_with_inds(X.column(i), temp);
            top_k(temp, R.column(i), I.column(i), TComp<indexed_value<T> >() );
        }
        
        return std::make_pair(mR, mI);     
    }
    else  // rows
    {
        marray mR = create_marray<T>(m, K);
        marray mI = create_marray<double>(m, K);
        array1d<indexed_value<T> > temp(n);
        
        const_aview2d<T, column_major_t> X = view2d<T>(mX);
        aview2d<T, column_major_t> R = view2d<T>(mR);
        aview2d<double, column_major_t> I = view2d<double>(mI);
        
        for (size_t i = 0; i < m; ++i)
        {
            export_with_inds(X.row(i), temp);
            top_k(temp, R.row(i), I.row(i), TComp<indexed_value<T> >() );
        }
        
        return std::make_pair(mR, mI);
    }             
}



template<typename T, template<typename U> class TComp>
inline marray do_top_k(const_marray mX, size_t K, int adim)
{    
    size_t m = mX.nrows();
    size_t n = mX.ncolumns();
        
    if (adim == 0)
    {        
        if (m == 1)
        {
            return vec_top_k<T, TComp>(mX, K, 2);
        }
        else if (n == 1)
        {
            return vec_top_k<T, TComp>(mX, K, 1);
        }
        else 
        {
            return mat_top_k<T, TComp>(mX, K, 1);
        }
    }
    else
    {   
        if (adim == 1 && n == 1)
        {
            return vec_top_k<T, TComp>(mX, K, 1);
        }
        else if (adim == 2 && m == 1)
        {
            return vec_top_k<T, TComp>(mX, K, 2);
        }
        else
        {
            return mat_top_k<T, TComp>(mX, K, adim);
        }
    }
}


template<typename T, template<typename U> class TComp>
inline std::pair<marray, marray> do_indexed_top_k(const_marray mX, size_t K, int adim)
{    
    size_t m = mX.nrows();
    size_t n = mX.ncolumns();
        
    if (adim == 0)
    {
        if (m == 1)
        {
            return vec_indexed_top_k<T, TComp>(mX, K, 2);
        }
        else if (n == 1)
        {
            return vec_indexed_top_k<T, TComp>(mX, K, 1);
        }
        else 
        {
            return mat_indexed_top_k<T, TComp>(mX, K, 1);
        }
    }
    else
    {   
        if (adim == 1 && n == 1)
        {
            return vec_indexed_top_k<T, TComp>(mX, K, 1);
        }
        else if (adim == 2 && m == 1)
        {
            return vec_indexed_top_k<T, TComp>(mX, K, 2);
        }        
        else
        {
            return mat_indexed_top_k<T, TComp>(mX, K, adim);
        }
    }
}



template<template<typename U> class TComp>
inline marray do_top_k_dtype(const_marray mX, size_t K, int adim)
{    
    switch (mX.class_id())
    {
        case mxDOUBLE_CLASS:
            return do_top_k<double, TComp>(mX, K, adim);
            
        case mxSINGLE_CLASS:
            return do_top_k<float, TComp>(mX, K, adim);
            
        case mxINT32_CLASS:
            return do_top_k<int32_t, TComp>(mX, K, adim);
            
        case mxUINT32_CLASS:
            return do_top_k<uint32_t, TComp>(mX, K, adim);
            
        default:
            throw mexception("top_k:invalidarg",
                "The X should be of type double, single, int32, or uint32.");
    }    
}

template<template<typename U> class TComp>
inline std::pair<marray, marray> do_indexed_top_k_dtype(const_marray mX, size_t K, int adim)
{    
    switch (mX.class_id())
    {
        case mxDOUBLE_CLASS:
            return do_indexed_top_k<double, TComp>(mX, K, adim);
            
        case mxSINGLE_CLASS:
            return do_indexed_top_k<float, TComp>(mX, K, adim);
            
        case mxINT32_CLASS:
            return do_indexed_top_k<int32_t, TComp>(mX, K, adim);
            
        case mxUINT32_CLASS:
            return do_indexed_top_k<uint32_t, TComp>(mX, K, adim);
            
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
 * [3]:  dim:   the dimension along which top K elements are selected
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
    const_marray mDim(prhs[3]);
    
    size_t K = (size_t)mK.get_scalar<double>();
    int code = (int)mCode.get_scalar<double>();
    int adim = (int)mDim.get_scalar<double>();
    
    marray mR, mI;
    
    if (nlhs <= 1)
    {
        if (code == 1)
        {
            mR = do_top_k_dtype<std::greater>(mX, K, adim);
        }
        else
        {
            mR = do_top_k_dtype<std::less>(mX, K, adim);
        }
        
        plhs[0] = mR.mx_ptr();
    }
    else
    {
        if (code == 1)
        {
            rbind(mR, mI) = do_indexed_top_k_dtype<std::greater>(mX, K, adim);
        }
        else 
        {
            rbind(mR, mI) = do_indexed_top_k_dtype<std::less>(mX, K, adim);
        }
        
        
        plhs[0] = mR.mx_ptr();
        plhs[1] = mI.mx_ptr();
    }
               
}


BCSMEX_MAINDEF
      
        
