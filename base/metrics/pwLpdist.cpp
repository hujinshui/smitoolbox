/********************************************************************
 *
 *  The C++ mex implementation of pairwise Minkowski distance,
 *  including L1-distance, L2-distance, and L-inf distance, as
 *  well as generic Minkowski distance.
 *
 *  Created by Dahua Lin, on Jun 2, 2008
 *
 ********************************************************************/


#include "mex.h"

#include <stdlib.h>
#include <math.h>

// type-dependent absolute value computation

template<typename Ty> struct calc_absdiff;

template<>
struct calc_absdiff<int> {
    inline int operator() (int x, int y) { return abs(x - y); }
};

template<>
struct calc_absdiff<float> {
    inline float operator() (float x, float y) { return fabs(x - y); }
};

template<>
struct calc_absdiff<double> {
    inline double operator() (double x, double y) { return fabs(x - y); }
};

template<typename Ty>
        inline Ty sqr(Ty x) {
    return x * x;
}



// p-dependent Computing routines

template<typename Ty>
void L1_dists(const void *pX1, const void *pX2, const void *pD, int n1, int n2, int d) {
    const Ty *X1 = (const Ty*)pX1;
    const Ty *X2 = (const Ty*)pX2;
    Ty *D = (Ty*)pD;
    
    const Ty *v1 = X1;
    const Ty *v2 = X2;
    
    calc_absdiff<Ty> calcabs;
    
    for (int j = 0; j < n2; ++j, v2 += d) {
        v1 = X1;
        for (int i = 0; i < n1; ++i, v1 += d) {
            Ty s = 0;
            
            for (int k = 0; k < d; ++k) {
                s += calcabs(v1[k], v2[k]);
            }
            
            *D++ = s;
        }
    }
}

template<typename Ty>
void L2_dists(const void *pX1, const void *pX2, const void *pD, int n1, int n2, int d) {
    const Ty *X1 = (const Ty*)pX1;
    const Ty *X2 = (const Ty*)pX2;
    Ty *D = (Ty*)pD;
    
    const Ty *v1 = X1;
    const Ty *v2 = X2;
    
    for (int j = 0; j < n2; ++j, v2 += d) {
        v1 = X1;
        for (int i = 0; i < n1; ++i, v1 += d) {
            Ty s = 0;
            
            for (int k = 0; k < d; ++k) {
                s += sqr(v1[k] - v2[k]);
            }
            
            *D++ = sqrt(s);
        }
    }
}

template<typename Ty>
void Linf_dists(const void *pX1, const void *pX2, const void *pD, int n1, int n2, int d) {
    const Ty *X1 = (const Ty*)pX1;
    const Ty *X2 = (const Ty*)pX2;
    Ty *D = (Ty*)pD;
    
    const Ty *v1 = X1;
    const Ty *v2 = X2;
    
    calc_absdiff<Ty> calcabs;
    
    for (int j = 0; j < n2; ++j, v2 += d) {
        v1 = X1;
        for (int i = 0; i < n1; ++i, v1 += d) {
            Ty s = calcabs(*v1, *v2);
            
            for (int k = 1; k < d; ++k) {
                Ty t = calcabs(v1[k], v2[k]);
                if (t > s) s = t;
            }
            
            *D++ = s;
        }
    }
}


template<typename Ty>
        void Lp_dists(const void *pX1, const void *pX2, const void *pD, int n1, int n2, int d, Ty p) {
    const Ty *X1 = (const Ty*)pX1;
    const Ty *X2 = (const Ty*)pX2;
    Ty *D = (Ty*)pD;
    
    const Ty *v1 = X1;
    const Ty *v2 = X2;
    
    calc_absdiff<Ty> calcabs;
    Ty invp = (Ty)1 / p;
    
    for (int j = 0; j < n2; ++j, v2 += d) {
        v1 = X1;
        for (int i = 0; i < n1; ++i, v1 += d) {
            Ty s = 0;
            
            for (int k = 0; k < d; ++k) {
                s += pow(calcabs(v1[k], v2[k]), p);
            }
            
            *D++ = pow(s, invp);
            
        }
    }
}



/*
 * entry
 *
 * Input
 *    - X1, X2: the input matrices d x n1, and d x n2
 *    - p:      the norm-index (double scalar)
 *
 * Output
 *    - dists:  n1 x n2 pairwise distance matrix
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // take input
    
    const mxArray *mxX1 = prhs[0];
    const mxArray *mxX2 = prhs[1];
    const mxArray *mxp = prhs[2];
    
    int n1 = mxGetN(mxX1);
    int n2 = mxGetN(mxX2);
    int d = mxGetM(mxX1);
    
    mxClassID cid = mxGetClassID(mxX1);
    
    // check validity
    
    if (mxIsSparse(mxX1) || mxIsSparse(mxX2) || mxIsComplex(mxX1) || mxIsComplex(mxX2)) {
        mexErrMsgIdAndTxt("pwLpdist:invalidarg", 
                "X1 and X2 should be both real full matrices.");
    }
    
    if (mxGetM(mxX2) != d) {
        mexErrMsgIdAndTxt("pwLpdist:invalidarg",
                "X1 and X2 should have the same number of rows.");
    }
    
    if (mxGetClassID(mxX2) != cid) {
        mexErrMsgIdAndTxt("pwLpdist:invalidarg", 
                "X1 and X2 should be of the same value type.");
    }
    
    double p = mxGetScalar(mxp);
    bool is_inf = mxIsInf(p);
    
    if (p < 1) {
        mexErrMsgTxt("p should be a positive value with p >= 1");
    }
    
    if (p == 1 || is_inf) {
        if (!(cid == mxDOUBLE_CLASS || cid == mxSINGLE_CLASS || cid == mxINT32_CLASS)) {
            mexErrMsgIdAndTxt("pwLpdist:invalidarg",
                    "When p is 1 or inf, only double, single, and int32 types are supported.");
        }
    }
    else {
        if (!(cid == mxDOUBLE_CLASS || cid == mxSINGLE_CLASS)) {
            mexErrMsgIdAndTxt("pwLpdist:invalidarg", 
                    "when p is a value other than 1 or inf, only double and single types are supported.");
        }
    }
    
    // delegate
    
    mxArray *mxD = 0;
    
    if (cid == mxDOUBLE_CLASS) {
        mxD = mxCreateDoubleMatrix(n1, n2, mxREAL);
    }
    else {
        mxD = mxCreateNumericMatrix(n1, n2, cid, mxREAL);
    }
    
    if (n1 > 0 && n2 > 0 && d > 0) {
        if (p == 1) {
            switch (cid) {
                case mxDOUBLE_CLASS:
                    L1_dists<double>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d);
                    break;
                case mxSINGLE_CLASS:
                    L1_dists<float>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d);
                    break;
                case mxINT32_CLASS:
                    L1_dists<int>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d);
                    break;
            }
        }
        else if (p == 2) {
            switch (cid) {
                case mxDOUBLE_CLASS:
                    L2_dists<double>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d);
                    break;
                case mxSINGLE_CLASS:
                    L2_dists<float>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d);
                    break;
            }
        }
        else if (is_inf) {
            switch (cid) {
                case mxDOUBLE_CLASS:
                    Linf_dists<double>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d);
                    break;
                case mxSINGLE_CLASS:
                    Linf_dists<float>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d);
                    break;
                case mxINT32_CLASS:
                    Linf_dists<int>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d);
                    break;
            }
        }
        else {
            switch (cid) {
                case mxDOUBLE_CLASS:
                    Lp_dists<double>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d, p);
                    break;
                case mxSINGLE_CLASS:
                    Lp_dists<float>(mxGetData(mxX1), mxGetData(mxX2), mxGetData(mxD), n1, n2, d, (float)p);
                    break;
            }
        }
    }
    
    // output
    plhs[0] = mxD;
    
}
