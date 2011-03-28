/********************************************************************
 *
 *  pieces_cimp.cpp
 *
 *  The C++ mex implementation for pieces
 *
 *  Created by Dahua Lin, on Mar 22, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

using namespace bcs;
using namespace bcs::matlab;

template<typename T>
void do_pieces_left_u(const T* edges, const T* x, double* y, int n, int m)
{    
    for (int i = 0; i < n; ++i)
    {
        int k = 0;
        T v = x[i];
        while (k <= m && v >= edges[k]) ++k;        
        y[i] = k;        
    }
}

template<typename T>
void do_pieces_right_u(const T* edges, const T* x, double* y, int n, int m)
{
    for (int i = 0; i < n; ++i)
    {
        int k = 0;
        T v = x[i];
        while (k <= m && v > edges[k]) ++k;
        y[i] = k;
    }
}


template<typename T>
void do_pieces_left_s(const T* edges, const T* x, double* y, int n, int m)
{
    int k = 0;
    for (int i = 0; i < n; ++i)
    {
        T v = x[i];
        while (k <= m && v >= edges[k]) ++k;        
        y[i] = k;
    }
}

template<typename T>
void do_pieces_right_s(const T* edges, const T* x, double* y, int n, int m)
{
    int k = 0;
    for (int i = 0; i < n; ++i)
    {
        T v = x[i];
        while (k <= m && v > edges[k]) ++k;        
        y[i] = k;
    }
}


template<typename T>
inline void do_pieces(const T* edges, const T* x, double *y, int n, int m, 
        bool is_left, bool is_sorted)
{
    if (is_left) 
    {
        if  (is_sorted)
            do_pieces_left_s<T>(edges, x, y, n, m);
        else
            do_pieces_left_u<T>(edges, x, y, n, m);
    }
    else 
    {
        if (is_sorted)
            do_pieces_right_s<T>(edges, x, y, n, m);
        else
            do_pieces_right_u<T>(edges, x, y, n, m);
    }
}



/**
 * The main entry
 *
 *  Input:
 *      [0] x:          the values [double or single matrix]
 *      [1] edges:      the edges of the bins [1 x (m+1) or (m+1) x 1]
 *      [2] dir:        whether it is left-close: 1-left close, 0-right close
 *      [3] sorted:     whether the input x is sorted
 *                      if true, it uses a more efficient implementation O(n + m)
 *                      if not, it uses a standard implementation O(n m) 
 *  Output:
 *      [0] y:          the output (int32, 1 x n)
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_marray mX(prhs[0]);
    const_marray mEdges(prhs[1]);
    const_marray mDir(prhs[2]);
    const_marray mSorted(prhs[3]);
        
    bool is_left = mDir.get_scalar<bool>();
    bool is_sorted = mSorted.get_scalar<bool>();
        
            
    // main
        
    int nr = (int)mX.nrows();
    int nc = (int)mX.ncolumns();
    int n = nr * nc;
    int m = (int)mEdges.nelems() - 1;
    
    marray mY = create_marray<double>(nr, nc);
    double *y = mY.data<double>();
        
    if (mX.is_double())
    {        
        do_pieces<double>(mEdges.data<double>(), mX.data<double>(), 
                y, n, m, is_left, is_sorted);
    }
    else if (mX.is_double())
    {                
        do_pieces<float>(mEdges.data<float>(), mX.data<float>(), 
                y, n, m, is_left, is_sorted);
    }    
    
    // output
    plhs[0] = mY.mx_ptr();
    
}


BCSMEX_MAINDEF
        



