/********************************************************************
 *
 *  qpps.cpp
 *
 *  The C++ mex implementation of qpps.m
 *
 *  Created by Dahua Lin, on Jul 21, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

#include <algorithm>

using namespace bcs;
using namespace bcs::matlab;


template <typename T>
struct entry
{
    int i;
    T v;    
        
    bool operator < (const entry& rhs) const
    {
        return v > rhs.v;   // just want to sort it in descending order
    }    
    
    void set(int _i, T _v)
    {
        i = _i;
        v = _v;
    }
};


template<typename T>
void qpps(caview1d<T> f, aview1d<T> x, entry<T> *es)
{
    int d = (int)f.nelems();
    
    // pre-condition: x has all zeros
    
    // make and sort entries
        
    for (int i = 0; i < d; ++i) 
    {
        es[i].set(i, f[i]);
    }
    
    std::sort(es, es + d);
    
    
    // main
    
    int k = 0;
    T s = T(0);
    T cumsum = T(0);
    
    while (k < d)                
    {
        cumsum += es[k].v;        
        T new_s = cumsum - (k+1) * es[k].v;
        
        if (new_s < 1)
        {
            s = new_s;
            ++k;
        }
        else break;
    }
    
    T a = (1 - s) / k - es[k - 1].v;
            
    // calculate output (x)
    
    for (int i = 0; i < k; ++i)
    {
        x[es[i].i] = es[i].v + a;            
    }
        
}





template<typename T>
inline marray do_qpps(const_marray mF)
{
    index_t m = mF.nrows();
    index_t n = mF.ncolumns();
    
    marray mX;
    entry<T> *es = 0;
    
    if (m == 1 || n == 1)
    {
        if (m == 1)
        {
            mX = create_marray<T>(1, (size_t)n);
            es = new entry<T>[n];
        }
        else
        {
            mX = create_marray<T>((size_t)m, 1);
            es = new entry<T>[m];
        }
        
        caview1d<T> f = view1d<T>(mF);
        aview1d<T> x = view1d<T>(mX);
        
        qpps(f, x, es);        
    }
    else
    {
        mX = create_marray<T>((size_t)m, (size_t)n);
        
        caview2d<T, column_major_t> F = view2d<T>(mF); 
        aview2d<T, column_major_t> X = view2d<T>(mX);
        
        es = new entry<T>[m];
        
        for (index_t i = 0; i < n; ++i)
        {
            qpps(F.column(i), X.column(i), es);
        }
    }
    
    delete[] es;

    return mX;
}




/**********
 *
 *  main entry:
 *
 *  input:
 *    [0]: f  [d x n]
 *
 *  output:
 *    [0]: x  [d x n]
 *
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    if (nrhs != 1)
    {
        throw mexception("qpps:invalidarg", 
            "The number of input arguments should be 1 for qpps");
    }
    
    const_marray mF(prhs[0]);
    
    if (!( mF.is_matrix() && mF.is_float() && !mF.is_sparse() && !mF.is_complex() ))
    {
        throw mexception("qpps:invalidarg", 
            "F should be a non-sparse real matrix.");
    }
    
    
    // main
    
    marray mX;
    
    if (mF.is_double())
    {
        mX = do_qpps<double>(mF);
    }
    else
    {
        mX = do_qpps<float>(mF);
    }
    
    
    // output
    plhs[0] = mX.mx_ptr();
}


BCSMEX_MAINDEF





