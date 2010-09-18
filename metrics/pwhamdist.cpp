/********************************************************************
 *
 *  C++ mex implementation of pairwise hamming distance computing
 *
 *  Created by Dahua Lin, on Jun 3, 2008
 *
 ********************************************************************/

#include "mex.h"

template<typename T>
void pwhamdist(const T *X1, const T *X2, double *D, int d, int n1, int n2)
{
    const T *v1 = X1;
    const T *v2 = X2;

    for (int j = 0; j < n2; ++j, v2 += d)
    {
        v1 = X1;

        for (int i = 0; i < n1; ++i, v1 += d)
        {
            int s = 0;
            for (int k = 0; k < d; ++k)
            {
                if (v1[k] != v2[k]) ++s;
            }

            *D++ = s;
        }	
    }
}


/*
 * entry
 *
 * Input
 *    X1, X2:  input matrices (d x n1 and d x n2 logical or numerical)
 * Output
 *    D:       pairwise distance matrix (n1 x n2 int32)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    const mxArray *mxX1 = prhs[0];
    const mxArray *mxX2 = prhs[1];
    
    if (mxIsSparse(mxX1) || mxIsSparse(mxX2) || 
            mxIsComplex(mxX1) ||mxIsComplex(mxX2) ||
            mxGetNumberOfDimensions(mxX1) != 2 ||
            mxGetNumberOfDimensions(mxX2) != 2)
    {
        mexErrMsgIdAndTxt("pwhamdist:invalidarg",
                "X1 and X2 should be both non-sparse non-complex matrices.");
    }

    int n1 = mxGetN(mxX1);
    int n2 = mxGetN(mxX2);
    int d = mxGetM(mxX1);

    if (mxGetM(mxX2) != d)
    {
        mexErrMsgIdAndTxt("pwhamdist:invalidarg",
                "X1 and X2 should have the same number of rows.");
    }

    // compute
    
    mxClassID cid = mxGetClassID(mxX1);
    if (cid != mxGetClassID(mxX2))
    {
        mexErrMsgIdAndTxt("pwhamdist:invalidarg",
                "X1 and X2 should have the same value class.");
    }
    
    
    mxArray *mxD = mxCreateDoubleMatrix(n1, n2, mxREAL);
    
    switch (cid)
    {
        case mxDOUBLE_CLASS:
            pwhamdist((const double*)mxGetData(mxX1), (const double*)mxGetData(mxX2), 
                mxGetPr(mxD), d, n1, n2);
            break;
        case mxSINGLE_CLASS:
            pwhamdist((const float*)mxGetData(mxX1), (const float*)mxGetData(mxX2), 
                mxGetPr(mxD), d, n1, n2);
            break;
        case mxLOGICAL_CLASS:
            pwhamdist((const bool*)mxGetData(mxX1), (const bool*)mxGetData(mxX2), 
                mxGetPr(mxD), d, n1, n2);
            break;
        case mxINT32_CLASS:
            pwhamdist((const int*)mxGetData(mxX1), (const int*)mxGetData(mxX2), 
                mxGetPr(mxD), d, n1, n2);
            break;
        case mxUINT32_CLASS:
            pwhamdist((const unsigned int*)mxGetData(mxX1), (const unsigned int*)mxGetData(mxX2), 
                mxGetPr(mxD), d, n1, n2);
            break;
        default:
            mexErrMsgIdAndTxt("pwhamdist:invalidarg", 
                    "Only double, single, logical, int32 and uint32 types are supported.");
    }        
    

    // output

    plhs[0] = mxD;
}


