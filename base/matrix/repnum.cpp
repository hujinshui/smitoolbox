/********************************************************************
 *
 *  repnum_cimp.cpp
 *
 *  The C++ mex implementation for repnum
 *
 *  Created by Dahua Lin, on June 6, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/array/array_eval.h>

using namespace bcs;
using namespace bcs::matlab;


template<typename T0>
array1d<size_t> get_nums(const_marray mNs)
{
    size_t n = mNs.nelems();
    array1d<size_t> ns(n);
    
    const T0 *in = mNs.data<T0>();
    
    for (index_t i = 0; i < (index_t)n; ++i)
    {
        ns(i) = (size_t)(in[i]);
    }
    
    return ns;
}


array1d<size_t> do_get_nums(const_marray mNs)
{
    if (!(mNs.is_vector() && !mNs.is_sparse()))
    {
        throw mexception("repnum:invalidarg", 
            "ns should be a non-sparse numeric vector.");
    }
    
    switch (mNs.class_id())
    {
        case mxDOUBLE_CLASS:
            return get_nums<double>(mNs);
        case mxSINGLE_CLASS:
            return get_nums<float>(mNs);
        case mxINT32_CLASS:
            return get_nums<int32_t>(mNs);
        case mxUINT32_CLASS:
            return get_nums<uint32_t>(mNs);
        case mxINT16_CLASS:
            return get_nums<int16_t>(mNs);
        case mxUINT16_CLASS:
            return get_nums<uint16_t>(mNs); 
        case mxINT8_CLASS:
            return get_nums<int8_t>(mNs);
        case mxUINT8_CLASS:
            return get_nums<uint8_t>(mNs);            
        default:
            throw mexception("repnum:invalidarg", 
                "ns should be a non-sparse numeric vector.");
    }
}


void check_vs_size(const_marray mNs, const_marray mVs)
{
    if (!(mVs.is_vector() && !mVs.is_sparse()))
    {
        throw mexception("repnum:invalidarg", 
            "vs should be a non-sparse numeric vector.");
    }
    
    if (mVs.nrows() != mNs.nrows() || mVs.ncolumns() != mVs.ncolumns())
    {
        throw mexception("repnum:invalidarg", 
            "The sizes of vs and ns are inconsistent.");
    }    
}


marray do_repnum(const_marray mNs)
{
    array1d<size_t> ns = do_get_nums(mNs);
    
    size_t n = ns.nelems();
    size_t N = sum(ns);
    
    marray mR = mNs.nrows() == 1 ? 
        create_marray<double>(1, N) :
        create_marray<double>(N, 1);
        
    double *r = mR.data<double>();
    
    for (size_t i = 0; i < n; ++i)
    {
        double v = i + 1;
        for (size_t j = 0; j < ns(i); ++j)
        {                        
            *(r++) = v;
        }
    }
    
    return mR;
}


template<typename T>
marray do_repnum_ex(const_marray mVs, const_marray mNs)
{
    check_vs_size(mNs, mVs);    
    array1d<size_t> ns = do_get_nums(mNs);
    const_aview1d<T> vs = view1d<T>(mVs);
    
    size_t n = ns.nelems();
    size_t N = sum(ns);
    
    marray mR = mNs.nrows() == 1 ? 
        create_marray<T>(1, N) :
        create_marray<T>(N, 1);
        
    T *r = mR.data<T>();
    
    for (size_t i = 0; i < n; ++i)
    {
        double v = vs(i);
        for (size_t j = 0; j < ns(i); ++j)
        {                        
            *(r++) = v;
        }
    }
    
    return mR;
    
}




/**
 *
 * main entry:
 *
 *  Input:
 *    [0] - ns 
 *
 *    or 
 *    [0] - vs, [1] - ns
 *
 *    ns:  the numbers of repeating times
 *    vs:  the values to be repeated
 *
 *  Output:
 *    [1]: r:   the generated vector 
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    marray mR;
    
    if (nrhs == 1)
    {
        const_marray mNs(prhs[0]);
        
        mR = do_repnum(mNs);
    }
    else if (nrhs == 2)
    {
        const_marray mVs(prhs[0]);
        const_marray mNs(prhs[1]);
        
        mR;
        
        switch (mVs.class_id())
        {
            case mxDOUBLE_CLASS:
                mR = do_repnum_ex<double>(mVs, mNs);
                break;
            case mxSINGLE_CLASS:
                mR = do_repnum_ex<float>(mVs, mNs);
                break;
            case mxINT32_CLASS:
                mR = do_repnum_ex<int32_t>(mVs, mNs);
                break;
            case mxUINT32_CLASS:
                mR = do_repnum_ex<uint32_t>(mVs, mNs);
                break;
            case mxLOGICAL_CLASS:
                mR = do_repnum_ex<bool>(mVs, mNs);
                break;
            case mxINT8_CLASS:
                mR = do_repnum_ex<int8_t>(mVs, mNs);
                break;
            case mxUINT8_CLASS:
                mR = do_repnum_ex<uint8_t>(mVs, mNs);
                break;
            case mxINT16_CLASS:
                mR = do_repnum_ex<int16_t>(mVs, mNs);
                break;
            case mxUINT16_CLASS:
                mR = do_repnum_ex<uint16_t>(mVs, mNs);
                break;
            default:
                throw mexception("repnum:invalidarg",
                    "vs should be of numeric or logical type");                
        }
    }
    else
    {
        throw mexception("repnum:invalidarg", 
            "The number of inputs to repnum should be 1 or 2.");
    }
    
    plhs[0] = mR.mx_ptr();
        
}

BCSMEX_MAINDEF



