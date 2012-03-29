/********************************************************************
 *
 *  aggreg_percol_cimp.cpp
 *
 *  The C++ mex implementation of aggreg.m
 *
 *  Created by Dahua Lin, on Feb 26, 2012
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


template<typename T, class Aggregator>
marray aggreg(const_marray mX, const_marray mI, int K)
{
    typedef typename Aggregator::result_type Tout;
    
    int m = mX.nrows();
    int n = mX.ncolumns();
    marray mR = create_marray<Tout>(K, n);
    
    Aggregator ag;
    
    const T *X = mX.data<T>();
    const int32_t *I = mI.data<int32_t>();    
    Tout *R = mR.data<Tout>();
    
    Tout v0 = ag.init();
    
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < K; ++i) R[i] = v0;
        
        for (int i = 0; i < m; ++i)
        {
            int32_t k = I[i];            
            if (k >= 0 && k < K)
            {
                ag(R[k], X[i]);
            }
        }
        
        
        X += m;
        I += m;
        R += K;
    }
    
    
    return mR;
}


template<typename T>
inline marray do_aggreg(const_marray mX, const_marray mI, int K, int code)
{        
    if (code == 1)
    {
        return aggreg<T, sum_ag<T> >(mX, mI, K);
    }
    else if (code == 2)
    {
        return aggreg<T, min_ag<T> >(mX, mI, K);
    }
    else // code == 3
    {
        return aggreg<T, max_ag<T> >(mX, mI, K);
    }
}


/**
 * main entry:
 *
 * Inputs:
 *  [0]:  X:    data [double|single|int32 matrix]
 *  [1]:  K:    #classes [double scalar]
 *  [2]:  I:    indices [zero-based int32]
 *  [3]:  code: function code (1-sum, 2-min, 3-max)
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
    const_marray mCode(prhs[3]);
     
    int K = (int)mK.get_scalar<double>();
    int code = (int)mCode.get_scalar<double>();
    
    // main delegate
    
    marray mR;
    
    if (mX.is_double())
    {
        mR = do_aggreg<double>(mX, mI, K, code);
    }
    else if (mX.is_single())
    {
        mR = do_aggreg<float>(mX, mI, K, code);
    }
    else if (mX.is_int32())
    {
        mR = do_aggreg<int32_t>(mX, mI, K, code);
    }
    else if (mX.is_logical())
    {
        mR = do_aggreg<bool>(mX, mI, K, code);        
    }
    else
    {
        throw mexception("aggreg_percol:invalidarg", 
            "aggreg only supports matrices of type double, single, int32, or logical");
    }
   
    // output
    plhs[0] = mR.mx_ptr();
    
}


BCSMEX_MAINDEF






