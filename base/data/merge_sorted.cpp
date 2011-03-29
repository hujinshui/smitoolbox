/********************************************************************
 *
 *  merge_sorted.cpp
 *
 *  The C++ mex implementation of merge_sorted
 *
 *  Created by Dahua Lin, on Mar 28, 2011
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

using namespace bcs;
using namespace bcs::matlab;


template<typename T>
marray do_merge_1(const mxArray *srcs[], int dim)
{
    const_marray mS(srcs[0]);    
    return duplicate(mS);
}


template<typename T>
marray do_merge_2(const mxArray *srcs[], int dim)
{
    const_marray mS1(srcs[0]);
    const_marray mS2(srcs[1]);
    
    const T *v1 = mS1.data<T>();
    const T *v2 = mS2.data<T>();
    
    size_t n1 = mS1.nelems();
    size_t n2 = mS2.nelems();
    
    marray mC = dim == 1 ? 
        create_marray<T>(n1+n2, 1) :
        create_marray<T>(1, n1+n2);
        
    T *c = mC.data<T>();
    
    size_t i1 = 0;
    size_t i2 = 0;
    
    while (i1 < n1 && i2 < n2)
    {
        if (v1[i1] <= v2[i2])
        {
            *c++ = v1[i1++];
        }
        else
        {
            *c++ = v2[i2++];
        }
    }
    
    while (i1 < n1)
    {
        *c++ = v1[i1++];
    }
    
    while (i2 < n2)
    {
        *c++ = v2[i2++];
    }    
        
    return mC;
}


template<typename T>
marray do_merge_multi(int K, const mxArray *srcs[], int dim)
{
    // extract input
    
    const T** ptrs = new const T*[K];
    size_t *ns = new size_t[K];    
    int *ps = new int[K];
    size_t N = 0;
    
    for (int k = 0; k < K; ++k)
    {
        const_marray mS(srcs[k]);
        
        ptrs[k] = mS.data<T>();
        ns[k] = mS.nelems();
        ps[k] = 0;
        
        N += ns[k];
    }
    
    
    // prepare output
    
    marray mC = dim == 1 ?
        create_marray<T>(N, 1) : 
        create_marray<T>(1, N);
        
    T *c = mC.data<T>();    
    
    // main loop
    
    for (size_t i = 0; i < N; ++i)
    {
        int sk = -1;
        T cv;
        
        for (int k = 0; k < K; ++k)
        {
            if (ps[k] < ns[k])
            {
                if (sk < 0 || ptrs[k][ps[k]] < cv)
                {
                    sk = k;
                    cv = ptrs[k][ps[k]];
                }                
            }
        }
        
        *(c++) = cv;
        ++ ps[sk];                
    }
    
    // release
    
    delete[] ptrs;
    delete[] ns;
    delete[] ps;
    
    return mC;
}



template<typename T>
marray do_merge_sorted(int K, const mxArray *srcs[], int dim)        
{
    if (K == 1)
    {
        return do_merge_1<T>(srcs, dim);
    }
    else if (K == 2)
    {
        return do_merge_2<T>(srcs, dim);
    }
    else
    {
        return do_merge_multi<T>(K, srcs, dim);
    }    
}


/**
 * main entry:
 *
 * Input:  sorted vectors
 *
 * Output: combined vectors
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    marray mR;
    
    if (nrhs > 0)
    {
        int K = nrhs;
        
        int dim = 0;
        mxClassID cid;
        
        // check inputs and determine dim
        
        for (int k = 0; k < K; ++k)
        {
            const_marray mS(prhs[k]);
            
            if (!( mS.is_vector() && mS.is_numeric() && !mS.is_sparse() ))
            {
                throw mexception("merge_sorted:invalidarg", 
                    "each input should be a numeric vector.");
            }
            
            if (k == 0)
            {
                cid = mS.class_id();                                
            }
            else
            {
                if (mS.class_id() != cid)
                {
                    throw mexception("merge_sorted:invalidarg", 
                        "The types of input vectors are inconsistent.");
                }
            }
            
            if (dim == 0)
            {
                if (mS.nrows() > 1)
                {
                    dim = 1;
                }
                else if (mS.ncolumns() > 1)
                {
                    dim = 2;
                }
            }                       
        }
        
        if (dim == 0) dim = 2; // by default
        
        // main delegate                
        
        switch (cid)
        {
            case mxDOUBLE_CLASS:
                mR = do_merge_sorted<double>(K, prhs, dim);
                break;
                
            case mxSINGLE_CLASS:
                mR = do_merge_sorted<float>(K, prhs, dim);
                break;                
                
            case mxINT32_CLASS:
                mR = do_merge_sorted<int32_t>(K, prhs, dim);
                break;        
                
            case mxUINT32_CLASS:
                mR = do_merge_sorted<uint32_t>(K, prhs, dim);
                break;
                
            case mxINT16_CLASS:
                mR = do_merge_sorted<int16_t>(K, prhs, dim);
                break;        
                
            case mxUINT16_CLASS:
                mR = do_merge_sorted<uint16_t>(K, prhs, dim);
                break; 
                
            case mxINT8_CLASS:
                mR = do_merge_sorted<int8_t>(K, prhs, dim);
                break;        
                
            case mxUINT8_CLASS:
                mR = do_merge_sorted<uint8_t>(K, prhs, dim);
                break;      
                
            default:
                throw mexception("merge_sorted:invalidarg", 
                    "Unsupported numeric class encountered.");
        }
        
    }
    else
    {
        mR = create_marray<double>((size_t)0, (size_t)0);
    }
    
    plhs[0] = mR.mx_ptr();
}


BCSMEX_MAINDEF




