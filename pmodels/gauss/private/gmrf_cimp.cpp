/********************************************************************
 *
 *  gmrf_cimp.cpp
 *
 *  A C++ mex to support some operations related to gmrf
 *
 *  Created by Dahua Lin, on Feb 12, 2011
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

using namespace bcs;
using namespace bcs::matlab;


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
    
    size_t m = (size_t)(const_marray(prhs[1]).get_scalar<double>());
    const_marray mS(prhs[2]);
    const_marray mX(prhs[3]);
    
    const int32_t *s = mS.data<int32_t>();
    
    size_t n = mX.nrows();
    size_t K = mX.ncolumns();
    
    marray mV = create_marray<double>(m, K);
    
    if (K == 1)
    {
        const_aview1d<double> x = view1d<double>(mX);
        aview1d<double> v = view1d<double>(mV);
        
        for (index_t i = 0; i < (index_t)m; ++i)
        {
            v(i) = x(s[i]);
        }
    }
    else
    {
        const_aview2d<double, column_major_t> X = view2d<double>(mX);
        aview2d<double, column_major_t> V = view2d<double>(mV); 
        
        for (index_t k = 0; k < K; ++k)
        {                       
            for (index_t i = 0; i < m; ++i)
            {
                V(i, k) = X(s[i], k);
            }
        }
    }
    
    plhs[0] = mV.mx_ptr();
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
    
    size_t m = (size_t)(const_marray(prhs[1]).get_scalar<double>());
    const_marray mS(prhs[2]);
    const_marray mT(prhs[3]);
    const_marray mX(prhs[4]);
    
    const int32_t *s = mS.data<int32_t>();
    const int32_t *t = mT.data<int32_t>();
    
    size_t n = mX.nrows();
    size_t K = mX.ncolumns();
    
    marray mV = create_marray<double>(m, K);
    
    if (K == 1)
    {
        const_aview1d<double> x = view1d<double>(mX);
        aview1d<double> v = view1d<double>(mV);
        
        for (index_t i = 0; i < m; ++i)
        {
            v(i) = TOp::calc(x(s[i]), x(t[i]));
        }
    }
    else
    {
        const_aview2d<double, column_major_t> X = view2d<double>(mX);
        aview2d<double, column_major_t> V = view2d<double>(mV); 
        
        for (index_t k = 0; k < K; ++k)
        {                      
            const_aview1d<double> x = X.column(k);
            
            for (index_t i = 0; i < m; ++i)
            {
                V(i, k) = TOp::calc(x(s[i]), x(t[i]));
            }
        }
    }
    
    plhs[0] = mV.mx_ptr();
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
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int code = (int)(const_marray(prhs[0]).get_scalar<double>());
    
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


BCSMEX_MAINDEF


