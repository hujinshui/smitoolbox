// Some special function codes

#ifndef SMI_SPECFUNCS_H
#define SMI_SPECFUNCS_H

#include <cmath>

namespace smi
{
    
    double gammaln(double x)
    {
        const double M_lnSqrt2PI = 0.91893853320467274178;
        
        const double s0 = 76.18009172947146;        
        const double s1 = -86.50532032941677;
        const double s2 = 24.01409824083091;
        const double s3 = -1.231739572450155;
        const double s4 = 0.1208650973866179e-2;
        const double s5 = -0.5395239384953e-5;
        
        /* Lanczos method */
        double denom = x+1;
        double x1 = x + 5.5;
        
        double series = 1.000000000190015 + 
                s0 / (x+1) + 
                s1 / (x+2) + 
                s2 / (x+3) + 
                s3 / (x+4) +
                s4 / (x+5) + 
                s5 / (x+6);
        
        return( M_lnSqrt2PI + 
                (x+0.5)* std::log(x1) - x1 + std::log(series/x) );
    }
    
       
    double digamma(double x)
    {
        const double c = 12;
        const double d1 = -0.57721566490153286;
        const double d2 = 1.6449340668482264365; /* pi^2/6 */
        const double s = 1e-6;
        const double s3 = 1./12;
        const double s4 = 1./120;
        const double s5 = 1./252;
        const double s6 = 1./240;
        const double s7 = 1./132;
        
        /* Use Taylor series if argument <= S */
        if(x <= s)
        {
            return d1 - 1/x + d2*x;
        }
        
        /* Reduce to digamma(X + N) where (X + N) >= c */
        double y = 0;
        while(x < c)
        {
            y -= 1.0 / x;
            x += 1.0;
        }
        
        /* Use de Moivre's expansion if argument >= C */
        /* This expansion can be computed in Maple via asympt(Psi(x),x) */
        if(x >= c)
        {
            double r = 1.0 / x;
            y += log(x) - 0.5*r;
            r *= r;
            
            double t = (s5 - r * (s6 - r * s7));
            y -= r * (s3 - r * (s4 - r * t));
        }
        return y;
    }


}

#endif
