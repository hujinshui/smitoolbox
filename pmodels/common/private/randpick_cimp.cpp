/********************************************************************
 *
 *  C++ mex implementation of core algorithm of randpick
 *
 *  Created by Dahua Lin, on Aug 10, 2011
 *
 *******************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <set>
#include <cmath>

using namespace bcs;
using namespace bcs::matlab;

mxArray *call_rand(int n, int len)
{
    mxArray *x = 0;
    
    mxArray *prhs[2];
    prhs[0] = mxCreateDoubleScalar(1);
    prhs[1] = mxCreateDoubleScalar(len);
    
    mexCallMATLAB(1, &x, 2, prhs, "rand");
    
    mxDestroyArray(prhs[0]);
    mxDestroyArray(prhs[1]);
    
    return x;
}

mxArray *call_randi(int n, int len)
{
    mxArray *x = 0;
    
    mxArray *prhs[3];
    prhs[0] = mxCreateDoubleScalar(n);
    prhs[1] = mxCreateDoubleScalar(1);
    prhs[2] = mxCreateDoubleScalar(len);
    
    mexCallMATLAB(1, &x, 3, prhs, "randi");
    
    mxDestroyArray(prhs[0]);
    mxDestroyArray(prhs[1]);
    mxDestroyArray(prhs[2]);
    
    return x;
}

void randpick_by_rejection_sampling(int n, int k, double *r)
{        
    std::set<int> s;
    
    int len = k + (k / 2);
    int remain = k;
    
    int it = k;
    
    while(remain > 0 && it > 0)
    {            
        mxArray *mxX = call_randi(n, len);
        const double *x = mxGetPr(mxX);
        
        for(int i = 0; i < len; ++i)
        {
            int v = (int)x[i];
            
            if (s.find(v) == s.end()) // not in set
            {
                s.insert(v);
                *(r++) = v;
                
                if (--remain == 0) break;
            }            
        }
        
        -- it;
        
        mxDestroyArray(mxX);
    }
}

void randpick_by_random_suffling(int n, int k, double *r)
{    
    // initialize
        
    mxArray *mxX = call_rand(n, k);
    const double *x = mxGetPr(mxX);
    
    int *s = new int[n];
    for (int i = 0; i < n; ++i) s[i] = i;
       
    for (int i = 0; i < k; ++i)
    {
        // pick
        int j = i + std::floor(x[i] * (n - i));
        
        // swap
        if (i != j)
        {
            int t = s[i];
            s[i] = s[j];
            s[j] = t;
        }
    }
    
    // copy
    for (int i = 0; i < k; ++i)
    {
        r[i] = (double)(s[i] + 1);
    }
    
    delete[] s;
    
    mxDestroyArray(mxX);
}



/**
 * main entry:
 *
 * Input
 *   [0]: n:  the total number of population (double)
 *   [1]: k:  the number of samples to pick (double)
 *   [2]: b:  whether to use shuffle
 *
 * Ouput
 *   [0]: r:  the resultant sample vector (k x 1)
 */
void bcsmex_main(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const_marray mN(prhs[0]);
    const_marray mK(prhs[1]);
    const_marray mB(prhs[2]);
    
    int n = mN.get_scalar<int32_t>();
    int k = mK.get_scalar<int32_t>();
    bool b = mB.get_scalar<bool>();
    
    marray mR = create_marray<double>((size_t)k, 1);    
    double *r = mR.data<double>();
    
    if (b)
    {
        randpick_by_random_suffling(n, k, r);        
    }
    else
    {
        randpick_by_rejection_sampling(n, k, r);
    }
    
    plhs[0] = mR.mx_ptr();
}


BCSMEX_MAINDEF

        