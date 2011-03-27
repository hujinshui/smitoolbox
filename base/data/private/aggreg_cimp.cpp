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


template<class Aggregator, typename TIn, class TIndexerIn, typename TOut, class TIndexerOut>
inline void aggreg_vector(
        int n, int K, Aggregator ag, 
        const_aview1d<TIn, TIndexerIn> x, 
        const const_aview1d<int32_t>& I, 
        aview1d<TOut, TIndexerOut> r)
{
    for (int i = 0; i < n; ++i)
    {
        int k = I[i];
        if (k >= 0 && k < K)
        {
            ag(r[k], x[i]);
        }
    }
}

// do aggregation

template<typename T, class Aggregator>
marray do_aggreg_rows(const_marray mX, const_marray mI, size_t K, Aggregator ag)
{    
    index_t m = (index_t)mX.nrows();
    index_t n = (index_t)mX.ncolumns();
 
    typedef typename Aggregator::result_type TOut;    
    marray mR = create_marray<TOut>((size_t)m, K); 
    
    const_aview1d<int32_t> I = view1d<int32_t>(mI); 
    
    if (m == 1)
    {
        const_aview1d<T> x = view1d<T>(mX);
        aview1d<TOut> r = view1d<TOut>(mR);
        fill(r, ag.init());
        
        aggreg_vector(n, (int)K, ag, x, I, r);
    }
    else
    {
        const_aview2d<T, column_major_t> X = view2d<T>(mX);
        aview2d<TOut, column_major_t> R = view2d<TOut>(mR);
        fill(R, ag.init());
        
        for (index_t i = 0; i < m; ++i)
        {
            aggreg_vector(n, (int)K, ag, X.row(i), I, R.row(i));
        }
    }
   
    return mR;
}


template<typename T, class Aggregator>
marray do_aggreg_columns(const_marray mX, const_marray mI, size_t K, Aggregator ag)
{
    size_t m = mX.nrows();
    size_t n = mX.ncolumns();
    
    typedef typename Aggregator::result_type TOut;    
    marray mR = create_marray<TOut>((size_t)K, n); 
    
    const_aview1d<int32_t> I = view1d<int32_t>(mI); 
    
    if (n == 1)
    {
        const_aview1d<T> x = view1d<T>(mX);
        aview1d<TOut> r = view1d<TOut>(mR);
        fill(r, ag.init());
        
        aggreg_vector(m, (int)K, ag, x, I, r);
    }
    else
    {
        const_aview2d<T, column_major_t> X = view2d<T>(mX);
        aview2d<TOut, column_major_t> R = view2d<TOut>(mR);
        fill(R, ag.init());
        
        for (index_t i = 0; i < n; ++i)
        {
            aggreg_vector(m, (int)K, ag, X.column(i), I, R.column(i));
        }
    }
   
    return mR;
}


template<typename T>
inline marray do_aggreg(const_marray mX, const_marray mI, size_t K, int code)
{    
    if (mI.nrows() != 1)
    {
        switch (code)
        {
            case SUM_CODE:
                return do_aggreg_columns<T, sum_ag<T> >(mX, mI, K, sum_ag<T>());
            case MIN_CODE:
                return do_aggreg_columns<T, min_ag<T> >(mX, mI, K, min_ag<T>());
            case MAX_CODE:
                return do_aggreg_columns<T, max_ag<T> >(mX, mI, K, max_ag<T>());                
        }        
    }
    else
    {
        switch (code)
        {
            case SUM_CODE:
                return do_aggreg_rows<T, sum_ag<T> >(mX, mI, K, sum_ag<T>());
            case MIN_CODE:
                return do_aggreg_rows<T, min_ag<T> >(mX, mI, K, min_ag<T>());
            case MAX_CODE:
                return do_aggreg_rows<T, max_ag<T> >(mX, mI, K, max_ag<T>());                
        }
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
     
    size_t K = (size_t)mK.get_scalar<double>();
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
        throw mexception("aggreg:invalidarg", 
            "aggreg only supports matrices of type double, single, int32, or logical");
    }
   
    // output
    plhs[0] = mR.mx_ptr();
    
}


BCSMEX_MAINDEF
        
        

