/********************************************************************
 *
 *  valueseg_cimp.cpp
 *
 *  The C++ mex implementation of valueseg
 *
 *  Created by Dahua Lin, on May 26, 2010
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <vector>

using namespace bcs;
using namespace bcs::matlab;


template<typename T>
void find_seg(size_t n, const T *v, std::vector<int>& offsets)
{
    offsets.push_back(0);    
    for (size_t i = 1; i < n; ++i)
    {
        if (v[i] != v[i-1]) offsets.push_back(i);
    }    
}


void do_make_segs(int n, const std::vector<int>& offsets, bool make_row, 
        marray& mSp, marray& mEp)
{        
    int m = (int)offsets.size() - 1;
        
    mSp = make_row ? create_marray<double>(1, m+1) : create_marray<double>(m+1, 1);        
    mEp = make_row ? create_marray<double>(1, m+1) : create_marray<double>(m+1, 1);
        
    double *sp = mSp.data<double>();
    double *ep = mEp.data<double>();
            
    for (int i = 0; i < m; ++i)
    {
        sp[i] = offsets[i] + 1;
        ep[i] = offsets[i + 1];
    }
    
    sp[m] = offsets[m] + 1;
    ep[m] = n;   
}


/***********
 *
 * Main entry:
 * 
 * Input:
 *   [0]: v:  the value vector to segment (numeric or logical or char)
 * Output:
 *   [0]: sp: the starting positions 
 *   [1]: ep: the ending positions
 *
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_marray mV(prhs[0]);
    
    if (!(mV.is_vector() && !mV.is_sparse() && !mV.is_empty()))
    {
        throw mexception("valueseg:invalidarg", 
            "x should be a full non-empty vector.");
    }
    
    size_t n = mV.nelems();
    bool make_row = mV.nrows() == 1;
        
    
    // delegate by type
        
    mxClassID cid = mV.class_id();
    
    std::vector<int> offsets;
    
    switch (cid)
    {
        case mxDOUBLE_CLASS:
            find_seg(n, mV.data<double>(), offsets);
            break;
        case mxSINGLE_CLASS:
            find_seg(n, mV.data<float>(), offsets);
            break;
        case mxINT32_CLASS:
            find_seg(n, mV.data<int32_t>(), offsets);
            break;
        case mxCHAR_CLASS:
            find_seg(n, mV.data<char>(), offsets);
            break;
        case mxUINT8_CLASS:
            find_seg(n, mV.data<uint8_t>(), offsets);
            break;
        case mxUINT32_CLASS:
            find_seg(n, mV.data<uint32_t>(), offsets);
            break;
        case mxLOGICAL_CLASS:
            find_seg(n, mV.data<bool>(), offsets);
            break;
        case mxINT8_CLASS:
            find_seg(n, mV.data<int8_t>(), offsets);
            break;
        case mxINT16_CLASS:
            find_seg(n, mV.data<int16_t>(), offsets);
            break;
        case mxUINT16_CLASS:
            find_seg(n, mV.data<uint16_t>(), offsets);
            break;
        default:
            throw mexception("valueseg:invalidarg", 
                "valueseg only supports numeric, char, or logical types.");
    }
    
    marray mSp, mEp;
    do_make_segs(n, offsets, make_row, mSp, mEp);
    
    // output
    
    plhs[0] = mSp.mx_ptr();
    plhs[1] = mEp.mx_ptr();
    
}


BCSMEX_MAINDEF





