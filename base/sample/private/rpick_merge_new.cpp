/********************************************************************
 *
 *  rpick_merge_new.cpp
 *
 *  A core-part of incremental sampling without replacement for 
 *  merging new samples into existing ones.
 *
 *  Created by Dahua Lin, on Sep 27, 2010
 *
 ********************************************************************/

#include <mex.h>


// merge two sorted sequence into a sorted one
// (x2 will be remapped to a disjoint domain)
void merge_seq(int n1, const double *x1, int n2, const double *x2, double *y)
{
    int i1 = 0;
    int i2 = 0;
    int k = 0;
    
    while (i1 < n1 && i2 < n2)
    {
        int cv2 = x2[i2] + i1;
        
        if (x1[i1] <= cv2)
        {
            y[k++] = x1[i1++];
        }
        else
        {
            y[k++] = cv2;
            ++i2;
        }
    }
    
    while (i1 < n1)
    {
        y[k++] = x1[i1++];
    }
    
    while (i2 < n2)
    {
        y[k++] = x2[i2++] + i1;
    }
}



// main entry:
// Input:
//   [0]:  x1 (existing sorted sequence)
//   [1]:  x2 (new sorted sequence)
// Output:
//   [0]:   merged sequence
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxX1 = prhs[0];
    const mxArray *mxX2 = prhs[1];
    
    int n1 = (int)mxGetNumberOfElements(mxX1);
    int n2 = (int)mxGetNumberOfElements(mxX2);
    
    const double *x1 = mxGetPr(mxX1);
    const double *x2 = mxGetPr(mxX2);
    
    // prepare output
    
    mxArray *mxY = mxCreateDoubleMatrix(1, n1+n2, mxREAL);
    double *y = mxGetPr(mxY);
    
    // merge
    
    merge_seq(n1, x1, n2, x2, y);
    
    // output
    plhs[0] = mxY;
}


