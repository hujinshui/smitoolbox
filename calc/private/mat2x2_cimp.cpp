/*************************************************************************
 *
 *  mat2x2_cimp.cpp
 *
 *  The C++ mex implementation of routines for fast real matrix 
 *  computation for 2 x 2 matrices.
 *
 *  Implemented computation routines and their codes
 *  - det (01 for gm, 02 for sm)
 *  - inv (03 for gm, 04 for sm)
 *  - chol (05 only for sm)
 *  - eigs (06 only for sm)
 *  - sqrtm (07 only for sm)
 *  
 *  Created by Dahua Lin, on Apr 7, 2010
 *
 *************************************************************************/

#include <mex.h>
#include <math.h>

template<typename T> inline T sqr(T v) { return v * v; }

template<typename T> 
inline T det_gm(const T *a)
{
    return a[0] * a[3] - a[1] * a[2];
}

template<typename T>
void do_det_gm(const T *a, T *r, int n)
{
    if (n == 1)
    {
        *r = det_gm(a);
    }
    else
    {
        for (int i = 0; i < n; ++i, a += 4, ++r) *r = det_gm(a);
    }
}


template<typename T>
inline T det_sm(const T *a)
{
    return a[0] * a[2] - sqr(a[1]);
}

template<typename T>
void do_det_sm(const T *a, T *r, int n)
{
    if (n == 1)
    {
        *r = det_sm(a);
    }
    else
    {
        for (int i = 0; i < n; ++i, a += 3, ++r) *r = det_sm(a);
    }
}


template<typename T>
inline void inv_gm(const T *a, T *r)
{
    T k = T(1) / det_gm(a);
    
    r[0] = k * a[3];
    r[1] = k * (- a[1]);
    r[2] = k * (- a[2]);
    r[3] = k * a[0];
}

template<typename T>
void do_inv_gm(const T *a, T *r, int n)
{
    if (n == 1)
    {
        inv_gm(a, r);
    }
    else
    {
        for (int i = 0; i < n; ++i, a += 4, r += 4) inv_gm(a, r);
    }
}


template<typename T>
inline void inv_sm(const T *a, T *r)
{
    T k = T(1) / det_sm(a);
    
    r[0] = k * a[2];
    r[1] = k * (-a[1]);
    r[2] = k * a[0];
}

template<typename T>
void do_inv_sm(const T *a, T *r, int n)
{
    if (n == 1)
    {
        inv_sm(a, r);
    }
    else
    {
        for (int i = 0; i < n; ++i, a += 3, r += 3) inv_sm(a, r);
    }
}



template<typename T>
inline void chol_sm(const T *a, T *r)
{
    T x = sqrt(a[0]);
    T y = a[1] / x;
    T z = sqrt(a[2] - sqr(y));
    
    r[0] = x;
    r[1] = y;
    r[2] = z;
}

template<typename T>
void do_chol_sm(const T *a, T *r, int n)
{
    if (n == 1)
    {
        chol_sm(a, r);
    }
    else
    {
        for (int i = 0; i < n; ++i, a += 3, r += 3) chol_sm(a, r);
    }
}



/*******
 *
 *  This function produces a decomposition of the eigen-system
 *
 *  A = R' * D * R
 *
 *  where D = diag([p, q]);
 *        R = diag([cos(t) -sin(t); sin(t) cos(t)]);
 *
 *  The range of t is (-pi/2, pi/2]
 ******/

template<typename T>
inline void eigs_sm(const T *a, T* r)
{
    T a0 = a[0];
    T a1 = a[1];
    T a2 = a[2];
        
    if (a1 == 0)
    {
        if (a0 >= a2)
        {
            r[0] = a0;
            r[1] = a2;
            r[2] = 0;
        }
        else
        {
            r[0] = a2;
            r[1] = a0;
            r[2] = M_PI_2;
        }
    }
    else
    {
        T u, t;
        if (a0 != a2)
        {
            T y = -2 * a1;
            T x = a0 - a2;
            
            t = atan2(y, x);
            u = fabs(x) > fabs(y) ? (x / cos(t)) : (y / sin(t));
            
            t /= 2;
            u /= 2;
        }
        else
        {
            if (a1 > 0)
            {
                t = - M_PI_4;
                u = a1;
            }
            else
            {
                t = M_PI_4;
                u = -a1;
            }
        }
        
        T m = (a0 + a2) / 2;
        r[0] = m + u;
        r[1] = m - u;
        r[2] = t;
    }
}

template<typename T>
void do_eigs_sm(const T *a, T *r, int n)
{
    if (n == 1)
    {
        eigs_sm(a, r);
    }
    else
    {
        for (int i = 0; i < n; ++i, a += 3, r += 3) eigs_sm(a, r);
    }
}



template<typename T>
inline void sqrtm_sm(const T *a, T* r)
{    
    T a0 = a[0];
    T a1 = a[1];
    T a2 = a[2];
    
    if (a1 == 0)
    {
        r[0] = sqrt(a0);
        r[1] = 0;
        r[2] = sqrt(a2);
    }
    else if (a0 == a2)
    {
        T u = a0 + a2;
        T v = 2 * a1;
        T s = sqrt((u + v) / 2);
        T t = sqrt((u - v) / 2);        
        
        r[0] = r[2] = (s + t) / 2;
        r[1] = (s - t) / 2;
    }
    else
    {
        T y = -2 * a1;
        T x = a0 - a2;
            
        T t = atan2(y, x);               
        T u = (fabs(x) > fabs(y) ? (x / cos(t)) : (y / sin(t))) / 2;
            
        t /= 2;
        T c = cos(t);
        T s = sin(t);
        
        T cc = c * c;
        T ss = s * s;
        T cs = c * s;
        
        T m = (a0 + a2) / 2;
        T p = sqrt(m + u);
        T q = sqrt(m - u);
        
        r[0] = p * cc + q * ss;
        r[1] = (q - p) * cs;
        r[2] = q * cc + p * ss;
    }        
}

template<typename T>
void do_sqrtm_sm(const T *a, T *r, int n)
{
    if (n == 1)
    {
        sqrtm_sm(a, r);
    }
    else
    {
        for (int i = 0; i < n; ++i, a += 3, r += 3) sqrtm_sm(a, r);
    }
}


/*****
 *
 * Main entry:
 *
 * Input
 *  [0]: source matrix
 *  [1]: code
 * Output
 *  [0]: result matrix
 *
 *******/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxA = prhs[0];
    const mxArray *mxCode = prhs[1];
    
    int code = (int)mxGetScalar(mxCode);
    
    int nelems = mxGetNumberOfElements(mxA);    
    int q = (code == 1 || code == 3) ? 4 : 3;    
    int d = (code < 3) ? 1 : q;
    int n = nelems / q;
    
    // main
    
    mxArray *mxR = 0;
    
    if (mxIsDouble(mxA))
    {
        const double *A = (const double*)mxGetData(mxA);                
        
        mxR = mxCreateNumericMatrix(d, n, mxDOUBLE_CLASS, mxREAL);
        double *R = (double*)mxGetData(mxR);
        
        switch (code)
        {            
            case 1:
                do_det_gm(A, R, n);
                break;                
            case 2:
                do_det_sm(A, R, n);
                break;                
            case 3:
                do_inv_gm(A, R, n);                
                break;                
            case 4:
                do_inv_sm(A, R, n);
                break;                
            case 5:
                do_chol_sm(A, R, n);
                break;                
            case 6:
                do_eigs_sm(A, R, n);                
                break;                
            case 7:
                do_sqrtm_sm(A, R, n);
                break;
        }        
    }
    else  // is single
    {
        const float *A = (const float*)mxGetData(mxA);                
        
        mxR = mxCreateNumericMatrix(d, n, mxSINGLE_CLASS, mxREAL);
        float *R = (float*)mxGetData(mxR);
        
        switch (code)
        {            
            case 1:
                do_det_gm(A, R, n);
                break;                
            case 2:
                do_det_sm(A, R, n);
                break;                
            case 3:
                do_inv_gm(A, R, n);                
                break;                
            case 4:
                do_inv_sm(A, R, n);
                break;                
            case 5:
                do_chol_sm(A, R, n);
                break;                
            case 6:
                do_eigs_sm(A, R, n);                
                break;                
            case 7:
                do_sqrtm_sm(A, R, n);
                break;
        }
    }
    
    // output
    
    plhs[0] = mxR;    
}











