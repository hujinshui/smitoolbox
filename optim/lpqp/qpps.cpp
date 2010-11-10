/********************************************************************
 *
 *  qpps.cpp
 *
 *  The C++ mex implementation of qpps.m
 *
 *  Created by Dahua Lin, on Jul 21, 2010
 *
 ********************************************************************/

#include <mex.h>

#include <algorithm>

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
void qpps(int d, const T *f, T *x, entry<T> *es)
{
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
void qpps_all(int d, int n, const T *f, T *x)
{
    // prepare work space
    
    entry<T> *es = new entry<T>[d];
    
    // main loop
    
    for (int i = 0; i < n; ++i, f += d, x += d)
    {
        qpps(d, f, x, es);
    }
    
    
    // release memory
    
    delete[] es;
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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    if (nrhs != 1)
        mexErrMsgIdAndTxt("qpps:invalidarg", "The number of input arguments should be 1 for qpps");
    
    const mxArray *mxF = prhs[0];
    
    if (!(mxGetNumberOfDimensions(mxF) == 2 && !mxIsSparse(mxF) && !mxIsComplex(mxF)))
        mexErrMsgIdAndTxt("qpps:invalidarg", "f must be a non-sparse real matrix");
    
    // main
    
    int d = mxGetM(mxF);
    int n = mxGetN(mxF);
    mxArray *mxX = 0;
    
    if (mxIsDouble(mxF))
    {
        const double *f = (const double*)mxGetData(mxF);
        mxX = mxCreateNumericMatrix(d, n, mxDOUBLE_CLASS, mxREAL);
        double *x = (double*)mxGetData(mxX);
        
        qpps_all(d, n, f, x);        
    }
    else if (mxIsSingle(mxF))
    {
        const float *f = (const float*)mxGetData(mxF);
        mxX = mxCreateNumericMatrix(d, n, mxSINGLE_CLASS, mxREAL);
        float *x = (float*)mxGetData(mxX);
        
        qpps_all(d, n, f, x);
    }
    else
    {
        mexErrMsgIdAndTxt("qpps:invalidarg", "f must be of double or single class.");
    }
    
    // output
    plhs[0] = mxX;
}








