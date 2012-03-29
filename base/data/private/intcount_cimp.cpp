/********************************************************************
 *
 *  intcount_cimp.cpp
 *
 *  The C++ mex implementation of integer counting
 *
 *  Created by Dahua Lin, on May 26, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

using namespace bcs;
using namespace bcs::matlab;


template<typename T>
inline void count_vec(int v0, int v1, const T *data, int len, double *c)
{
    for (size_t i = 0; i < len; ++i)
    {
        int cv = (int)(data[i]);
        if (cv >= v0 && cv <= v1)
        {
            ++ c[cv - v0];
        }
    }
}


template<typename T>
marray count(int v0, int v1, const_marray mA, bool b)
{
    int K = v1 - v0 + 1;    
    const T *data = !mA.is_empty() ? mA.data<T>() : NULL;
    
    if (!b)
    {
        marray mC = create_marray<double>(1, K);
        double *c = mC.data<double>();
        count_vec(v0, v1, data, mA.nelems(), c);
        return mC;
    }
    else
    {
             
        int m = (int)mA.nrows();
        int n = (int)mA.ncolumns();
        
        marray mC = create_marray<double>(K, n);
        double *c = mC.data<double>();
        
        for (int j = 0; j < n; ++j, c += K, data += m)
        {
            count_vec(v0, v1, data, m, c);
        }
        
        return mC;        
    }    
}


inline marray do_count(int v0, int v1, const_marray mA, bool b)
{            
    switch (mA.class_id())
    {
        case mxDOUBLE_CLASS:
            return count<double>(v0, v1, mA, b);
            
        case mxSINGLE_CLASS:
            return count<float>(v0, v1, mA, b);
            
        case mxINT32_CLASS:
            return count<int32_t>(v0, v1, mA, b);
            
        case mxUINT32_CLASS:
            return count<uint32_t>(v0, v1, mA, b);
            
        case mxINT16_CLASS:
            return count<int16_t>(v0, v1, mA, b);
            
        case mxUINT16_CLASS:
            return count<uint16_t>(v0, v1, mA, b);     
            
        case mxINT8_CLASS:
            return count<int8_t>(v0, v1, mA, b);
            
        case mxUINT8_CLASS:
            return count<uint8_t>(v0, v1, mA, b);
            
        default:                    
            throw mexception("intcount:invalidarg", 
                "Unsupported value type for V.");
    }
    
}


/***********
 *
 *  Main entry
 *   
 *  Input:
 *    [0]: v0:      the lower-end of the integer range [int32]
 *    [1]: v1:      the higher-end of the integer range [int32]
 *    [2]: A:       the value matrix/array [numeric]
 *    [3]: b:       whether to do the counting per-column [logical]
 *
 *  Output:
 *    [0]: C:       the matrix of counts
 *
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check take inputs  
    
    const_marray mV0(prhs[0]);
    const_marray mV1(prhs[1]);
    const_marray mA(prhs[2]);
    const_marray mB(prhs[3]);
    
    int v0 = mV0.get_scalar<int32_t>();
    int v1 = mV1.get_scalar<int32_t>();
    bool b = mB.get_scalar<bool>();
    
    
    // main
        
    plhs[0] = do_count(v0, v1, mA, b).mx_ptr();

}


BCSMEX_MAINDEF


