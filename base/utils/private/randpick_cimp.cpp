/********************************************************************
 *
 *  C++ mex implementation of core algorithm of randpick
 *
 *  Created by Dahua Lin, on Aug 10, 2011
 *
 *******************************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include <random>
#include <set>

using namespace bcs;
using namespace bcs::matlab;

typedef std::mt19937 rand_engine_t;

inline int randi(rand_engine_t& eng, int m)
{
    return (int)(eng() % (unsigned)m);
}


void randpick_by_rejection_sampling(int n, int k, double *r, rand_engine_t& eng)
{
    std::set<int> s;
    
    for (int i = 0; i < k; ++i)
    {
        int x = randi(eng, n);
        
        while (s.find(x) != s.end())
        {
            x = randi(eng, n);
        }
        
        s.insert(x);
        r[i] = (double)(x+1);
    }
}

void randpick_by_random_suffling(int n, int k, double *r, rand_engine_t& eng)
{
    // initialize
    int *s = new int[n];
    for (int i = 0; i < n; ++i) s[i] = i;
        
    for (int i = 0; i < k; ++i)
    {
        // pick
        int j = i + randi(eng, n-i);
        
        // swap
        int t = s[i];
        s[i] = s[j];
        s[j] = t;
    }
    
    // copy
    for (int i = 0; i < k; ++i)
    {
        r[i] = (double)(s[i] + 1);
    }
    
    delete[] s;
}



/**
 * main entry:
 *
 * Input
 *   [0]: n:  the total number of population (double)
 *   [1]: k:  the number of samples to pick (double)
 *   [2]: a:  algorithm to use (double)
 *            0 - rejection sampling
 *            1 - random suffling 
 *   [3]: seed: a random seed
 *
 * Ouput
 *   [0]: r:  the resultant sample vector (k x 1)
 */
void bcsmex_main(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const_marray mxN(prhs[0]);
    const_marray mxK(prhs[1]);
    const_marray mxA(prhs[2]);
    const_marray mxSeed(prhs[3]);
    
    int n = (int)mxN.get_scalar<double>();
    int k = (int)mxK.get_scalar<double>();
    int a = (int)mxA.get_scalar<double>();
    unsigned long seed = (unsigned long)mxSeed.get_scalar<uint32_t>();
    
    marray mxR = create_marray<double>((size_t)k, 1);    
    double *r = mxR.data<double>();
    
    rand_engine_t eng;
    eng.seed(seed);
    
    if (a == 0)
    {
        randpick_by_rejection_sampling(n, k, r, eng);
    }
    else
    {
        randpick_by_random_suffling(n, k, r, eng);
    }
    
    plhs[0] = mxR.mx_ptr();
}


BCSMEX_MAINDEF

        