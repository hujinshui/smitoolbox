/********************************************************************
 *
 *  valueseg_cimp.cpp
 *
 *  The C++ mex implementation of valueseg
 *
 *  Created by Dahua Lin, on May 26, 2010
 *
 ********************************************************************/

#include <mex.h>

#include <vector>

template<typename T>
void find_seg(int n, const T *v, std::vector<int>& offsets)
{
    offsets.push_back(0);    
    for (int i = 1; i < n; ++i)
    {
        if (v[i] != v[i-1]) offsets.push_back(i);
    }    
}


void mxMakeSegs(int n, const std::vector<int>& offsets, bool make_row, 
        mxArray*& mxSp, mxArray*& mxEp)
{        
    int m = offsets.size() - 1;
    
    mxSp = make_row ? mxCreateDoubleMatrix(1, m+1, mxREAL) :
        mxCreateDoubleMatrix(m+1, 1, mxREAL);
        
    mxEp = make_row ? mxCreateDoubleMatrix(1, m+1, mxREAL) :
        mxCreateDoubleMatrix(m+1, 1, mxREAL);
        
    double *sp = mxGetPr(mxSp);
    double *ep = mxGetPr(mxEp);
            
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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxV = prhs[0];
    
    int n = mxGetNumberOfElements(mxV);
    bool make_row = mxGetM(mxV) == 1;
    
    // delegate by type
    
    mxArray *mxSp = 0;
    mxArray *mxEp = 0;
    
    mxClassID cid = mxGetClassID(mxV);
    
    std::vector<int> offsets;
    const char *v = (char*)mxGetData(mxV);
    
    switch (cid)
    {
        case mxDOUBLE_CLASS:
            find_seg(n, (const double*)(v), offsets);
            break;
        case mxSINGLE_CLASS:
            find_seg(n, (const float*)(v), offsets);
            break;
        case mxINT32_CLASS:
            find_seg(n, (const int*)(v), offsets);
            break;
        case mxCHAR_CLASS:
            find_seg(n, (const char*)(v), offsets);
            break;
        case mxUINT8_CLASS:
            find_seg(n, (const char*)(v), offsets);
            break;
        case mxUINT32_CLASS:
            find_seg(n, (const unsigned int*)(v), offsets);
            break;
        case mxLOGICAL_CLASS:
            find_seg(n, (const bool*)(v), offsets);
            break;
        case mxINT8_CLASS:
            find_seg(n, (const unsigned char*)(v), offsets);
            break;
        case mxINT16_CLASS:
            find_seg(n, (const short*)(v), offsets);
            break;
        case mxUINT16_CLASS:
            find_seg(n, (const unsigned short*)(v), offsets);
            break;
        default:
            mexErrMsgIdAndTxt("valueseg:invalidarg", 
                    "valueseg only supports numeric, char, or logical types.");
    }
            
    mxMakeSegs(n, offsets, make_row, mxSp, mxEp);
    
    // output
    
    plhs[0] = mxSp;
    plhs[1] = mxEp;
    
}






