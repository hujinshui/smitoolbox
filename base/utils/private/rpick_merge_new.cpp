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

#include <bcslib/matlab/bcs_mex.h>

using namespace bcs;
using namespace bcs::matlab;


// merge two sorted sequence into a sorted one
// (x2 will be remapped to a disjoint domain)
void merge_seq(size_t n1, const double *x1, size_t n2, const double *x2, double *y)
{
    size_t i1 = 0;
    size_t i2 = 0;
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
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const_marray mX1(prhs[0]);
    const_marray mX2(prhs[1]);
    
    size_t n1 = mX1.nelems();
    size_t n2 = mX2.nelems();
    
    const double *x1 = mX1.data<double>();
    const double *x2 = mX2.data<double>();
    
    // prepare output
    
    marray mY = create_marray<double>(1, n1+n2);
    double *y = mY.data<double>();
    
    // merge
    
    merge_seq(n1, x1, n2, x2, y);
    
    // output
    plhs[0] = mY.mx_ptr();
}


BCSMEX_MAINDEF

        
