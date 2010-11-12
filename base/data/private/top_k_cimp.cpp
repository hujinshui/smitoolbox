/********************************************************************
 *
 *  top_k_cimp.cpp
 *
 *  The C++ mex implementation of top_k.m
 *
 *  Created by Dahua Lin, on Nov 11, 2010
 *
 ********************************************************************/

#include <mex.h>
#include <string.h>
#include <algorithm>
#include <functional>


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


template<typename T, template<typename U> class TComp>
inline void top_k(const T *x, int n, int K, T *r, T *ws)
{
    if (n > 1)
    {
        TComp<T> comp;
        
        ::memcpy(ws, x, n * sizeof(T));    
        std::nth_element(ws, ws+(K-1), ws+n, comp);
        std::sort(ws, ws+K, comp);
        ::memcpy(r, ws, K * sizeof(T));
    }
    else
    {
        TComp<T> comp;
        
        *r = *(std::min_element(x, x+n, comp));
    }    
}


template<typename T, template<typename U> class TComp>
inline void indexed_top_k(const T *x, int n, int K, 
        T *r, double *inds, indexed_value<T> *ws)
{
    if (n > 1)
    {        
        for (int i = 0; i < n; ++i)
        {
            ws[i].set(i, x[i]);
        }
        
        TComp<indexed_value<T> > comp;
        
        std::nth_element(ws, ws+(K-1), ws+n, comp);
        std::sort(ws, ws+K, comp);
        
        for (int i = 0; i < K; ++i)
        {
            r[i] = ws[i].v;
            inds[i] = ws[i].i + 1;
        }        
    }
    else
    {
        TComp<T> comp;
        
        int i = std::min_element(x, x+n, comp) - x;
        *r = x[i];
        *inds = i+1;
    }
}


template<typename T, template<typename U> class TComp>
void vec_top_k(int n, int K, const T *x, T *r, double *inds)
{    
    if (inds == 0)
    {
        T *ws = new T[n];
        top_k<T, TComp>(x, n, K, r, ws);
        delete[] ws;
    }
    else 
    {
        indexed_value<T> *ws = new indexed_value<T>[n];
        indexed_top_k<T, TComp>(x, n, K, r, inds, ws);
        delete[] ws;
    }
    
}


template<typename T, template<typename U> class TComp>
void mat_top_k(int n, int nc, int K, const T *x, T *r, double *inds)
{    
    if (inds == 0)
    {
        T *ws = new T[n];
        for (int i = 0; i < nc; ++i) 
        {
            top_k<T, TComp>(x + i*n, n, K, r + i*K, ws);
        }
        delete[] ws;
    }
    else 
    {
        indexed_value<T> *ws = new indexed_value<T>[n];
        for (int i = 0; i < nc; ++i) 
        {
            indexed_top_k<T, TComp>(x + i*n, n, K, r + i*K, inds + i*K, ws);
        }
        delete[] ws;
    }
               
}



template<typename T, template<typename U> class TComp>
inline void do_top_k(const mxArray *mxX, int K, mxClassID cid,  
        int nlhs, mxArray *plhs[])
{
    int m = mxGetM(mxX);
    int n = mxGetN(mxX);
    const T *x = (const T*)mxGetData(mxX);
        
    mxArray *mxR = 0;
    mxArray *mxInds = 0;
    double *inds = 0;        
    
    if (m == 1)
    {
        mxR = mxCreateNumericMatrix(1, K, cid, mxREAL);
        T *r = (T*)mxGetData(mxR);        
        if (nlhs >= 2)
        {
            mxInds = mxCreateDoubleMatrix(1, K, mxREAL);
            inds = mxGetPr(mxInds);
        }        
        
        vec_top_k<T, TComp>(n, K, x, r, inds);
    }
    else if (n == 1)
    {
        mxR = mxCreateNumericMatrix(K, 1, cid, mxREAL);
        T *r = (T*)mxGetData(mxR);        
        if (nlhs >= 2)
        {
            mxInds = mxCreateDoubleMatrix(K, 1, mxREAL);
            inds = mxGetPr(mxInds);
        }        
        
        vec_top_k<T, TComp>(m, K, x, r, inds);
    }
    else  // m > 1 && n > 1
    {
        mxR = mxCreateNumericMatrix(K, n, cid, mxREAL);
        T *r = (T*)mxGetData(mxR);
        if (nlhs >= 2)
        {
            mxInds = mxCreateDoubleMatrix(K, n, mxREAL);
            inds = mxGetPr(mxInds);
        }
        
        mat_top_k<T, TComp>(m, n, K, x, r, inds);
    }
    
    plhs[0] = mxR;
    if (nlhs >= 2)
    {
        plhs[1] = mxInds;
    }    
}


template<template<typename U> class TComp>
inline void do_top_k_dtype(
        const mxArray *mxX, int K, int nlhs, mxArray *plhs[])
{
    mxClassID cid = mxGetClassID(mxX);
    
    switch (cid)
    {
        case mxDOUBLE_CLASS:
            do_top_k<double, TComp>(mxX, K, cid, nlhs, plhs);
            break;
        case mxSINGLE_CLASS:
            do_top_k<float, TComp>(mxX, K, cid, nlhs, plhs);
            break;
        case mxINT32_CLASS:
            do_top_k<int, TComp>(mxX, K, cid, nlhs, plhs);
            break;
        case mxUINT32_CLASS:
            do_top_k<unsigned int, TComp>(mxX, K, cid, nlhs, plhs);
            break;
        default:
            mexErrMsgIdAndTxt("top_k:invalidarg",
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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *mxX = prhs[0];
    const mxArray *mxK = prhs[1];
    const mxArray *mxCode = prhs[2];
    
    int K = (int)mxGetScalar(mxK);
    int code = (int)mxGetScalar(mxCode);
    
    if (code == 1)
    {
        do_top_k_dtype<std::greater>(mxX, K, nlhs, plhs);
    }
    else
    {
        do_top_k_dtype<std::less>(mxX, K, nlhs, plhs);
    }
            
}


