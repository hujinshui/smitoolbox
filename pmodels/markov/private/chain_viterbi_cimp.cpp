/**********************************************************
 *
 *  chain_viterbi_cimp.cpp
 *
 *  The C++ mex implementation of chain_viterbi
 *
 *  Created by Dahua Lin, on Feb 1, 2012
 *
 **********************************************************/


#include <bcslib/matlab/bcs_mex.h>
#include <cmath>

using namespace bcs;
using namespace bcs::matlab;


struct FullSO
{
    FullSO(int K, const double *B)
    : _K(K), _B(B)
    {   
    }
    
    double operator() (int i, int j) const
    {
        return _B[i + j * _K];
    }                
    
    int _K;
    const double *_B;
};


struct SimpleSO
{
    SimpleSO(const double *B)
    : b0(B[0]), b1(B[1]) 
    {
    }
            
    double operator() (int i, int j) const
    {
        return i == j ? b0 : b1;
    }
        
    double b0;
    double b1;
};




