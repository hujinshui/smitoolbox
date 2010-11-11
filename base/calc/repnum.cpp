/********************************************************************
 *
 *  repnum_cimp.cpp
 *
 *  The C++ mex implementation for repnum
 *
 *  Created by Dahua Lin, on June 6, 2010
 *
 ********************************************************************/

#include <mex.h>

template<typename T>
int get_nums(int n, const T *ns0, int *ns)
{
    int s = 0;
    for (int i = 0; i < n; ++i)
    {
        s += (ns[i] = (int)(ns0[i]));
    }
    return s;
}


int get_nums(const mxArray *mxNs, int n, int *ns)
{    
    mxClassID cid = mxGetClassID(mxNs);
    
    if (cid == mxDOUBLE_CLASS)
    {
        return get_nums<double>(n, (const double*)mxGetData(mxNs), ns);
    }
    else if (cid == mxSINGLE_CLASS)
    {
        return get_nums<float>(n, (const float*)mxGetData(mxNs), ns);
    }
    else if (cid == mxINT32_CLASS)
    {
        return get_nums<int>(n, (const int*)mxGetData(mxNs), ns);
    }
    else if (cid == mxUINT32_CLASS)
    {
        return get_nums<unsigned int>(n, (const unsigned int*)mxGetData(mxNs), ns);
    }
    else
    {
        mexErrMsgIdAndTxt("repnum:invalidarg", 
                "ns should be of type double, single, int32, or uint32");
    }            
}


void do_repnum(int n, const int *ns, double *r)
{
    for (int k = 0; k < n; ++k)
    {
        int c = ns[k];
        double v = k+1;
        
        for (int i = 0; i < c; ++i) *r++ = v;
    }
}


template<typename T>
void do_repvals(int n, const T *vs, const int *ns, T *r)
{
    for (int k = 0; k < n; ++k)
    {
        int c = ns[k];
        T v = vs[k];
        
        for (int i = 0; i < c; ++i) *r++ = v;
    }
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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *mxNs = 0;
    const mxArray *mxVs = 0;
    
    if (nrhs == 1)
    {
        mxNs = prhs[0];        
    }
    else if (nrhs == 2)
    {
        mxVs = prhs[0];
        mxNs = prhs[1];
    }
    else
    {
        mexErrMsgIdAndTxt("repnum:invalidarg", 
                "The number of inputs to repnum should be 1.");
    }

    // check ns
            
    int m_ns = mxGetM(mxNs);
    int n_ns = mxGetN(mxNs);   
    if (!(!mxIsSparse(mxNs) && !mxIsComplex(mxNs) && !mxIsEmpty(mxNs) && 
            mxGetNumberOfDimensions(mxNs) == 2 && (m_ns == 1 || n_ns == 1)))
        mexErrMsgIdAndTxt("repnum:invalidarg", 
                "ns should be a non-empty real full vector.");    
    int n = m_ns * n_ns;
    
    int *ns = new int[n];
    int N = get_nums(mxNs, n, ns);
        
    // check vs
    
    if (mxVs != 0)
    {
        int m_vs = mxGetM(mxVs);
        int n_vs = mxGetN(mxVs);
        if (!(!mxIsSparse(mxVs) && !mxIsComplex(mxVs) && 
                mxGetNumberOfDimensions(mxVs) == 2 && (m_ns == m_vs && n_ns == n_vs)))
            mexErrMsgIdAndTxt("repnum:invalidarg", 
                "vs should be a non-empty real full vector with the same size of ns.");
    }
    
    // main
    
    // determine the size of return
    int m_r, n_r;    
    if (m_ns == 1)
    {
        m_r = 1;
        n_r = N;
    }
    else
    {
        m_r = N;
        n_r = 1;
    }
        
    
    mxArray *mxR = 0;
    
    if (mxVs == 0)
    {
        mxR = mxCreateDoubleMatrix(m_r, n_r, mxREAL);
        do_repnum(n, ns, mxGetPr(mxR));
    }
    else
    {
        mxClassID cid = mxGetClassID(mxVs);
        switch (mxGetClassID(mxVs))
        {
            case mxDOUBLE_CLASS:
                mxR = mxCreateNumericMatrix(m_r, n_r, cid, mxREAL);
                do_repvals(n, (double*)mxGetData(mxVs), ns, (double*)mxGetData(mxR));
                break;
                
            case mxSINGLE_CLASS:
                mxR = mxCreateNumericMatrix(m_r, n_r, cid, mxREAL);
                do_repvals(n, (float*)mxGetData(mxVs), ns, (float*)mxGetData(mxR));
                break;  
                
            case mxINT32_CLASS:
                mxR = mxCreateNumericMatrix(m_r, n_r, cid, mxREAL);
                do_repvals(n, (int*)mxGetData(mxVs), ns, (int*)mxGetData(mxR));
                break;
                
            case mxUINT32_CLASS:
                mxR = mxCreateNumericMatrix(m_r, n_r, cid, mxREAL);
                do_repvals(n, (unsigned int*)mxGetData(mxVs), ns, (unsigned int*)mxGetData(mxR));
                break;
                
            case mxINT16_CLASS:
                mxR = mxCreateNumericMatrix(m_r, n_r, cid, mxREAL);
                do_repvals(n, (short*)mxGetData(mxVs), ns, (short*)mxGetData(mxR));
                break;
                
            case mxUINT16_CLASS:
                mxR = mxCreateNumericMatrix(m_r, n_r, cid, mxREAL);
                do_repvals(n, (unsigned short*)mxGetData(mxVs), ns, (unsigned short*)mxGetData(mxR));
                break;
                
            case mxINT8_CLASS:
                mxR = mxCreateNumericMatrix(m_r, n_r, cid, mxREAL);
                do_repvals(n, (char*)mxGetData(mxVs), ns, (char*)mxGetData(mxR));
                break;
                
            case mxUINT8_CLASS:
                mxR = mxCreateNumericMatrix(m_r, n_r, cid, mxREAL);
                do_repvals(n, (unsigned char*)mxGetData(mxVs), ns, (unsigned char*)mxGetData(mxR));
                break;
                
            case mxCHAR_CLASS:
                {
                    mwSize dims[2] = {m_r, n_r};
                    mxR = mxCreateCharArray(2, dims);
                    do_repvals(n, (mxChar*)mxGetData(mxVs), ns, (mxChar*)mxGetData(mxR));
                }
                break;
                
        }
    }

                    
    delete[] ns;
    
    // output
    
    plhs[0] = mxR;
        
}

