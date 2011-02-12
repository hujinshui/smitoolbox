/********************************************************************
 *
 *  gmrf_cimp.cpp
 *
 *  A C++ mex to support some operations related to gmrf
 *
 *  Created by Dahua Lin, on Feb 12, 2011
 *
 ********************************************************************/

#include "../../../graph/clib/graph_mex.h"

using namespace smi;


struct prod
{
    static double calc(double x, double y)
    {
        return x * y;
    }
};

struct dsqr
{
    static double calc(double x, double y)
    {
        double d = x - y;
        return d * d;
    }
};


/** 
 * Task: take values at edge end points
 *
 * Input:
 *      [1]:    m (double scalar)
 *      [2]:    i (the end point indices)
 *      [3]:    X [n x K]
 *
 * Ouput
 *      [0]:    V [m x K]
 *
 */
void do_take_values(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    int m = (int)mxGetScalar(prhs[1]);
    MArray mS(prhs[2]);
    MArray mX(prhs[3]);
    
    const int *s = mS.get_data<int>();
    const double *X = mX.get_data<double>();
    
    int n = mX.nrows();
    int K = mX.ncols();
    
    mxArray *mxV = mxCreateDoubleMatrix(m, K, mxREAL);
    double *V = mxGetPr(mxV);
    
    if (K == 1)
    {
        for (int i = 0; i < m; ++i)
        {
            V[i] = X[s[i]];
        }
    }
    else
    {
        for (int k = 0; k < K; ++k)
        {
            const double *x = X + n * k;
            double *v = V + m * k;           
            
            for (int i = 0; i < m; ++i)
            {
                v[i] = x[s[i]];
            }
        }
    }
    
    plhs[0] = mxV;
}



/**
 * Task: compute the values at edges 
 *
 * Input:
 *      [1]:    m (double scalar)
 *      [2]:    es [2m x 1 int32]
 *      [3]:    et [2m x 1 int32]
 *      [4]:    X [n x K]
 *
 * Output:
 *      [0]:    V [m x K]
 */
template<typename TOp>
void do_compute_at_edges(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    int m = (int)mxGetScalar(prhs[1]);
    MArray mS(prhs[2]);
    MArray mT(prhs[3]);
    MArray mX(prhs[4]);
    
    const int *s = mS.get_data<int>();
    const int *t = mT.get_data<int>();
    const double *X = mX.get_data<double>();
    
    int n = mX.nrows();
    int K = mX.ncols();
    
    mxArray *mxV = mxCreateDoubleMatrix(m, K, mxREAL);
    double *V = mxGetPr(mxV);
    
    if (K == 1)
    {
        for (int i = 0; i < m; ++i)
        {
            V[i] = TOp::calc(X[s[i]], X[t[i]]);
        }
    }
    else
    {
        for (int k = 0; k < K; ++k)
        {
            const double *x = X + n * k;
            double *v = V + m * k;           
            
            for (int i = 0; i < m; ++i)
            {
                v[i] = TOp::calc(x[s[i]], x[t[i]]);
            }
        }
    }
    
    plhs[0] = mxV;
}



/**
 * main entry
 * 
 * Input:
 *      [0] code:   0 - x(s) or x(e): do_take_values
 *                  1 - compute x(s) .* x(e): do_compute_at_edges<prod>
 *                  2 - compute (x(s) - x(e)).^2: do_compute_at_edges<dsqr>
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int code = (int)mxGetScalar(prhs[0]);
    
    if (code == 0)
    {
        do_take_values(nlhs, plhs, nrhs, prhs);
    }
    else if (code == 1)
    {
        do_compute_at_edges<prod>(nlhs, plhs, nrhs, prhs);
    }    
    else if (code == 2)
    {
        do_compute_at_edges<dsqr>(nlhs, plhs, nrhs, prhs);
    }
}

