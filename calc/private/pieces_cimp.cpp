/********************************************************************
 *
 *  pieces_cimp.cpp
 *
 *  The C++ mex implementation of pieces
 *
 *  Created by Dahua Lin, on Mar 22, 2010
 *
 ********************************************************************/

#include <mex.h>


template<typename T>
void do_pieces_left(const T* edges, const T* x, int* y, int n, int m)
{    
    int i = 0;
    int k = 0;
                     
    for(; k <= m; ++k)
    {
        T e = edges[k];        
        for(; x[i] < e && i < n; ++i) y[i] = k;
    }
     
    for (; i < n; ++i) y[i] = k;        
}

template<typename T>
void do_pieces_right(const T* edges, const T* x, int* y, int n, int m)
{
    int i = 0; 
    int k = 0;
    
    for (; k <= m; ++k)
    {
        T e = edges[k];
        for(; x[i] <= e && i < n; ++i) y[i] = k;
    }
    
    for (; i < n; ++i) y[i] = k;
}



/**
 * The main entry
 *
 *  Input:
 *      [0] x:          the values [1 x n or n x 1]
 *      [1] edges:      the edges of the bins [1 x (m+1) or (m+1) x 1]
 *      [2] is_left:    whether left-close or right-close [logical scalar]
 *  Output:
 *      [0] y:          the output (int32, 1 x n)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    const mxArray *mxX = prhs[0];
    const mxArray *mxEdges = prhs[1];
    const mxArray *mxIsLeft = prhs[2];
    
    int n = (int)mxGetNumberOfElements(mxX);
    int m = (int)mxGetNumberOfElements(mxEdges) - 1;
    bool is_left = mxIsLogicalScalarTrue(mxIsLeft);        
    
    // main
    mxArray *mxY = mxCreateNumericMatrix(1, n, mxINT32_CLASS, mxREAL);
    int *y = (int*)mxGetData(mxY);
    
    if (mxIsDouble(mxX))
    {
        const double *x = (const double*)mxGetData(mxX);
        const double *edges = (const double*)mxGetData(mxEdges);
        
        if (is_left) 
            do_pieces_left<double>(edges, x, y, n, m);
        else
            do_pieces_right<double>(edges, x, y, n, m);
    }
    else if (mxIsSingle(mxX))
    {
        const float *x = (const float*)mxGetData(mxX);
        const float *edges = (const float*)mxGetData(mxEdges);
        
        if (is_left)
            do_pieces_left<float>(edges, x, y, n, m);
        else
            do_pieces_right<float>(edges, x, y, n, m);
    }    
    
    // output
    plhs[0] = mxY;
    
}

