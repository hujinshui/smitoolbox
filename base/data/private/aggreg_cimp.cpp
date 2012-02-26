/********************************************************************
 *
 *  aggreg_cimp.cpp
 *
 *  The C++ mex implementation of aggreg.m
 *
 *  Created by Dahua Lin, on Nov 10, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <limits>

using namespace bcs;
using namespace bcs::matlab;

const int SUM_CODE = 1;
const int MIN_CODE = 2;
const int MAX_CODE = 3;


// aggregators

template<typename T>
struct sum_ag
{
    typedef T result_type;
    
    result_type init() const
    {
        return T(0);
    }
    
    void operator() (result_type& r, const T& x) const
    {
        r += x;
    }
};

template<>
struct sum_ag<bool>
{
    typedef double result_type;
    
    result_type init() const
    {
        return false;
    }
    
    void operator() (result_type& r, const bool& x) const
    {
        r += (double)x;
    }
};


template<typename T>
struct min_ag
{
    typedef T result_type;
    
    result_type init() const
    {
        if (std::numeric_limits<T>::has_infinity)
        {
            return std::numeric_limits<T>::infinity();
        }
        else
        {
            return std::numeric_limits<T>::max();
        }                
    }
    
    void operator() (result_type& r, const T& x) const
    {
        if (x < r) r = x;
    }
};

template<>
struct min_ag<bool>
{
    typedef bool result_type;
    
    result_type init() const
    {
        return true;
    }
    
    void operator() (result_type& r, const bool& x) const
    {
        r &= x;
    }
};



template<typename T>
struct max_ag
{
    typedef T result_type;
    
    result_type init() const
    {
        if (std::numeric_limits<T>::has_infinity)
        {
            return - std::numeric_limits<T>::infinity();
        }
        else
        {
            return std::numeric_limits<T>::min();
        }                
    }
    
    void operator() (result_type& r, const T& x) const
    {
        if (x > r) r = x;
    }
};


template<>
struct max_ag<bool>
{
    typedef bool result_type;
    
    result_type init() const
    {
        return false;
    }
    
    void operator() (result_type& r, const bool& x) const
    {
        r |= x;
    }
};


// core algorithm



template<typename T, class Aggregator>
marray aggreg(const_marray mX, const_marray mI, int dim, int K)
{
    typedef typename Aggregator::result_type Tout;
    
    Aggregator ag;
    Tout vinit = ag.init();
    
    int m = mX.nrows();
    int n = mX.ncolumns();
    
    const int32_t *I = mI.data<int32_t>();
    const T *x = mX.data<T>();
    
    
    if (dim == 1)
    {
        marray mR = create_marray<Tout>(K, n);
        Tout *r = mR.data<Tout>();
        
        int rn = K * n;
        for (int i = 0; i < rn; ++i) r[i] = vinit;                
    
        if (n == 1)
        {                        
            for (int i = 0; i < m; ++i)
            {
                int32_t k = I[i];
                if (k >= 0 && k < K) ag(r[k], x[i]);
            }
        }
        else
        {
            for (int j = 0; j < n; ++j, x += m, r += K)
            {
                for (int i = 0; i < m; ++i)
                {
                    int32_t k = I[i];
                    if (k >= 0 && k < K) ag(r[k], x[i]);
                }
            }
        }
        
        
        return mR;
    }
    else  // dim == 2
    {
        marray mR = create_marray<Tout>(m, K);
        Tout *R = mR.data<Tout>();
        
        int rn = m * K;
        for (int i = 0; i < rn; ++i) R[i] = vinit;
                
        if (m == 1)
        {
            for (int j = 0; j < n; ++j)
            {
                int32_t k = I[j];
                if (k >= 0 && k < K) ag(R[k], x[j]);
            }
        }
        else
        {
            for (int j = 0; j < n; ++j, x += m)
            {
                int32_t k = I[j];
                                
                if (k >= 0 && k < K)
                {
                    Tout *r = R + k * m;
                    
                    for (int i = 0; i < m; ++i) ag(r[i], x[i]);
                }
            }            
        }        
        
        return mR;
    }
    
}


template<typename T>
inline marray do_aggreg(const_marray mX, const_marray mI, int K, int dim, int code)
{    
    if (code == 1)
    {
        return aggreg<T, sum_ag<T> >(mX, mI, dim, K);
    }
    else if (code == 2)
    {
        return aggreg<T, min_ag<T> >(mX, mI, dim, K);
    }
    else // code == 3
    {
        return aggreg<T, max_ag<T> >(mX, mI, dim, K);
    }
}




/**
 * main entry:
 *
 * Inputs:
 *  [0]:  X:    data [double|single|int32 matrix]
 *  [1]:  K:    #classes [double scalar]
 *  [2]:  I:    indices [zero-based int32]
 *  [3]:  dim:  the direction along which the aggregation is performed
 *  [4]:  code: function code (1-sum, 2-min, 3-max)
 *
 * Outputs:
 *  [0]:  R:    the results
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_marray mX(prhs[0]);
    const_marray mK(prhs[1]);
    const_marray mI(prhs[2]);
    const_marray mDim(prhs[3]);
    const_marray mCode(prhs[4]);
     
    int K = (int)mK.get_scalar<double>();
    int dim = (int)mDim.get_scalar<double>();
    int code = (int)mCode.get_scalar<double>();
    
    // main delegate
    
    marray mR;
    
    if (mX.is_double())
    {
        mR = do_aggreg<double>(mX, mI, K, dim, code);
    }
    else if (mX.is_single())
    {
        mR = do_aggreg<float>(mX, mI, K, dim, code);
    }
    else if (mX.is_int32())
    {
        mR = do_aggreg<int32_t>(mX, mI, K, dim, code);
    }
    else if (mX.is_logical())
    {
        mR = do_aggreg<bool>(mX, mI, K, dim, code);        
    }
    else
    {
        throw mexception("aggreg:invalidarg", 
            "aggreg only supports matrices of type double, single, int32, or logical");
    }
   
    // output
    plhs[0] = mR.mx_ptr();
    
}


BCSMEX_MAINDEF
        
        

