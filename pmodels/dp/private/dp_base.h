/********************************************************************
 *
 *  dp_base.h
 *
 *  Basic facilities for DP-based implementation
 *
 *  Created by Dahua Lin, on Sep 17, 2011
 *
 ********************************************************************/

#ifndef DP_BASE_H_
#define DP_BASE_H_


/**
 * Draw from a discrete distribution 
 */
inline int dd_draw(int K, const double *w, double tw, double u)
{
    double v = u * tw;
    
    int k = -1;
    double cs = 0.0;
    int Km1 = K - 1;
        
    while (cs < v && k < Km1)
    {
        cs += w[++k];
    }
    
    if (cs < v) ++ k; // indicates the draw falls outside range
    
    return k;
}




#endif

