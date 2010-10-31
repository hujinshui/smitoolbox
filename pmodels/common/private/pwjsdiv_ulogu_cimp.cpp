/********************************************************************
 *
 *  C++ mex implementation of the u*log(u) part of 
 *  pairwise J-S divergence computation
 *
 *  History
 *      Created by Dahua Lin, on Jun 5, 2008
 *
 ********************************************************************/

#include "mex.h"

#include <math.h>

template<typename Ty>
void pwjsdiv_ulogu(const Ty *P, const Ty *Q, Ty *D, int d, int n1, int n2)
{
    const Ty *p = P;
    const Ty *q = Q;

    for (int j = 0; j < n2; ++j, q += d)
    {
	p = P;

	for (int i = 0; i < n1; ++i, p += d)
	{
	    Ty s = 0;

	    for (int k = 0; k < d; ++k)
	    {
		Ty pv = p[k];
		Ty qv = q[k];				
				
		if (pv > 0 || qv > 0)
		{		    
		    Ty u = (pv + qv) / 2;

		    s += u * log(u);
		} 
	    }

	    *D++ = s;
	}
    }
}

/*
 * entry
 *
 * Input
 *     P, Q:  the pmf matrices (d x n1 and d x n2)
 * Output
 *     D:     the pairwise divergence matrix (n1 x n2)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input

    const mxArray *mxP = prhs[0];
    const mxArray *mxQ = prhs[1];

    int d = mxGetM(mxP);
    int n1 = mxGetN(mxP);
    int n2 = mxGetN(mxQ);

    // check validity
    
    if (mxIsSparse(mxP) || mxIsComplex(mxP) || mxGetNumberOfDimensions(mxP) != 2)
    {
	mexErrMsgTxt("P should be a full real matrix.");
    }

    if (mxIsSparse(mxQ) || mxIsComplex(mxQ) || mxGetNumberOfDimensions(mxQ) != 2)
    {
	mexErrMsgTxt("Q should be a full real matrix.");
    }

    if (mxGetM(mxQ) != d)
    {
	mexErrMsgTxt("Q should have the same number of rows as in P.");
    }

    mxClassID cid = mxGetClassID(mxP);
    if (mxGetClassID(mxQ) != cid)
    {
	mexErrMsgTxt("P and Q should have the same value class.");
    }

    // main
    mxArray *mxD = 0;
    
    switch (cid)
    {
	case mxDOUBLE_CLASS:
	    mxD = mxCreateDoubleMatrix(n1, n2, mxREAL);
	    pwjsdiv_ulogu<double>(mxGetPr(mxP), mxGetPr(mxQ), mxGetPr(mxD), d, n1, n2);
	    break;

	case mxSINGLE_CLASS:
	    mxD = mxCreateNumericMatrix(n1, n2, mxSINGLE_CLASS, mxREAL);
	    pwjsdiv_ulogu<float>(
		    (const float*)mxGetData(mxP), 
		    (const float*)mxGetData(mxQ),
		    (float*)mxGetData(mxD), d, n1, n2);
	    break;

	default:
	    mexErrMsgTxt("Only double or single class is supported.");
    }

    // output

    plhs[0] = mxD;
}
