/********************************************************************
 *
 *  C++ mex implementation of X' * P with rule inf * 0 = 0
 *
 *  Created by Dahua Lin, on Jun 5, 2008
 *
 ********************************************************************/

#include "mex.h"

template<typename Ty>
void mmtimes(const Ty *P, const Ty *X, Ty *R, int d, int n1, int n2)
{
    const Ty *p = P;
    const Ty *x = X;

    Ty *ps = new Ty[d];
    int *nz = new int[d];

    for (int j = 0; j < n1; ++j, p += d)
    {
	// prepare for non-zero P's
	int c = 0;

	for (int k = 0; k < d; ++k)
	{
	    Ty pv = p[k];
	    if (pv > 0)
	    {
		ps[c] = pv;
		nz[c++] = k;
	    }
	}	
	
	// compute
	x = X;
	for (int i = 0; i < n2; ++i, x += d)
	{
	    Ty s = 0;

	    for (int k = 0; k < c; ++k)
	    {
		s += ps[k] * x[nz[k]];
	    }	    

	    *R++ = s;
	}
    }

    delete[] ps;
    delete[] nz;
}

/*
 * entry
 *
 * Input
 *     - P:  the measure (probability) d x n1
 *     - X:  the facility d x n2
 * Output
 *     - R:  the result, n2 x n1
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    const mxArray *mxP = prhs[0];
    const mxArray *mxX = prhs[1];

    int d = mxGetM(mxP);
    int n1 = mxGetN(mxP);
    int n2 = mxGetN(mxX);

    // check validity
    
    if (mxIsSparse(mxP) || mxIsComplex(mxP) || mxGetNumberOfDimensions(mxP) != 2)
    {
	mexErrMsgTxt("P should be a full real matrix.");
    }

    if (mxIsSparse(mxX) || mxIsComplex(mxX) || mxGetNumberOfDimensions(mxX) != 2)
    {
	mexErrMsgTxt("X should be a full real matrix.");
    }

    if (mxGetM(mxX) != d)
    {
	mexErrMsgTxt("X should have the same number of rows as in P.");
    }

    mxClassID cid = mxGetClassID(mxP);

    if (mxGetClassID(mxX) != cid)
    {
	mexErrMsgTxt("X and P should be of the same class.");
    }

    // main
    mxArray *mxD = 0;
    
    switch (cid)
    {
	case mxDOUBLE_CLASS:
	    mxD = mxCreateDoubleMatrix(n2, n1, mxREAL);
	    mmtimes(mxGetPr(mxP), mxGetPr(mxX), mxGetPr(mxD), d, n1, n2);
	    break;

	case mxSINGLE_CLASS:
	    mxD = mxCreateNumericMatrix(n2, n1, mxSINGLE_CLASS, mxREAL);
	    mmtimes(
		    (const float*)mxGetData(mxP), 
		    (const float*)mxGetData(mxX),
		    (float *)mxGetData(mxD), d, n1, n2);		    
	    break;

	default:
	    mexErrMsgTxt("P and X should be either double or single.");
    }

    // output
    plhs[0] = mxD;

}


