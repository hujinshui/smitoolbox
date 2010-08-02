/********************************************************************
 *
 *  C++ mex implementation of pairwise hamming distance computing
 *
 *  Created by Dahua Lin, on Jun 3, 2008
 *
 ********************************************************************/

#include "mex.h"

void pwhamdist(const bool *X1, const bool *X2, int *D, int d, int n1, int n2)
{
    const bool *v1 = X1;
    const bool *v2 = X2;

    for (int j = 0; j < n2; ++j, v2 += d)
    {
	v1 = X1;

	for (int i = 0; i < n1; ++i, v1 += d)
	{
	    int s = 0;
	    for (int k = 0; k < d; ++k)
	    {
		s += (int)(v1[k] ^ v2[k]);
	    }

	    *D++ = s;
	}	
    }
}


/*
 * entry
 *
 * Input
 *    X1, X2:  input matrices (d x n1 and d x n2 logical)
 * Output
 *    D:       pairwise distance matrix (n1 x n2 int32)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    const mxArray *mxX1 = prhs[0];
    const mxArray *mxX2 = prhs[1];

    int n1 = mxGetN(mxX1);
    int n2 = mxGetN(mxX2);
    int d = mxGetM(mxX1);

    // check validity

    if (mxIsSparse(mxX1) || mxIsSparse(mxX2))
    {
	mexErrMsgTxt("X1 and X2 should be both logical matrices.");
    }

    if (!mxIsLogical(mxX1) || !mxIsLogical(mxX2))
    {
	mexErrMsgTxt("X1 and X2 should be both logical matrices.");
    }

    if (mxGetM(mxX2) != d)
    {
	mexErrMsgTxt("X2 should have the same number of rows as in X1.");
    }

    // compute
    mxArray *mxD = mxCreateNumericMatrix(n1, n2, mxINT32_CLASS, mxREAL);
    
    pwhamdist((const bool*)mxGetData(mxX1), (const bool*)mxGetData(mxX2), 
	    (int*)mxGetData(mxD), d, n1, n2);

    // output

    plhs[0] = mxD;
}
