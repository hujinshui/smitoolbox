/********************************************************************
 *
 *  ddsample.h
 *
 *  The header files for sampling from discrete distribution
 *
 *  Created by Dahua Lin, on Nov 6, 2010
 *
 ********************************************************************/


#ifndef MSP_DDSAMPLE_H
#define MSP_DDSAMPLE_H

namespace msp
{

template<typename T>
inline T sum(int n, const T *x)
{
    T s(0);
    for (int i = 0; i < n; ++i) s += x[i];
    return s;
}   
    

template<typename T>
inline T cumsum(int n, const T *x, T *cs)
{
    T s(0);
    for (int i = 0; i < n; ++i)
        cs[i] = (s += x[i]);
    
    return s;
}


template<typename T>
inline int locate_interval(int n, const T *redges, T v)
{
    int i = 0;
    while (i < n && v > redges[i]) ++i;
    return i;
}


inline void count_labels(int K, int n, const int *labels, int *counts)
{
    for (int k = 0; k < K; ++k) counts[k] = 0;
    
    for (int i = 0; i < n; ++i)
    {
        int l = labels[i];
        if (l >= 0 && l < K) ++ counts[l];
    }
}


}

#endif

