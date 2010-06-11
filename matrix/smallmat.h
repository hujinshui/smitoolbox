/********************************************************************
 *
 *  smallmat.h
 *
 *  The common header shared by all small matrix computation routines
 *
 *  Created by Dahua Lin, on June 8, 2010
 *
 ********************************************************************/

#include <mex.h>
#include <cmath>

/// The struct to represent small matrices

template<typename T>
struct Vec2
{
    T v1;
    T v2;
};

template<typename T>
struct Vec3
{
    T v1;
    T v2;
    T v3;
};

template<typename T>
struct Mat2x2
{
    T v11;
    T v21;
    T v12;
    T v22;
};

template<typename T>
struct SMat2x2
{
    T v11;
    T v12;  
    T v22;
};


template<typename T>
struct Polar2
{
    T a;
    T b;
    T theta;
};

template<typename T>
struct LTMat2x2
{
    T v11;
    T v21;
    T v22;
};




// input manipulation

inline void general_input_check(const mxArray *mxA, const char *msgid)
{
    if (mxIsEmpty(mxA) || mxIsSparse(mxA) || mxIsComplex(mxA))
        mexErrMsgIdAndTxt(msgid,
            "The input array must be a non-empty non-sparse real array.");
}

inline bool test_size(const mxArray *mxA, int nr, int nc, int& n)
{
    if (nc == 1)
    {
        if (mxGetNumberOfDimensions(mxA) == 2 && mxGetM(mxA) == nr)
        {
            n = mxGetN(mxA);
            return true;
        }
        else
        {
            return false;
        }               
    }
    else
    {
        int ndim = mxGetNumberOfDimensions(mxA);
        if (ndim == 2)
        {
            if (mxGetM(mxA) == nr && mxGetN(mxA) == nc)
            {            
                n = 1;
                return true;
            }
            else
            {
                return false;
            }
        }
        else if (ndim == 3)
        {
            const mwSize *dims = mxGetDimensions(mxA);
            
            if (nr == dims[0] && nc == dims[1])
            {
                n = (int)(dims[2]);
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
}



template<typename T> mxArray *create_mat(int m, int n);
template<typename T> mxArray *create_cube(int m, int n, int p);


template<>
inline mxArray* create_mat<double>(int m, int n)
{
    return mxCreateDoubleMatrix(m, n, mxREAL);
}

template<>
inline mxArray* create_mat<float>(int m, int n)
{
    return mxCreateNumericMatrix(m, n, mxSINGLE_CLASS, mxREAL);
}

      
template<>
inline mxArray* create_cube<double>(int m, int n, int p)
{
    mwSize dims[3] = {(mwSize)m, (mwSize)n, (mwSize)p};
    return mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
}

template<>
inline mxArray* create_cube<float>(int m, int n, int p)
{
    mwSize dims[3] = {(mwSize)m, (mwSize)n, (mwSize)p};
    return mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
}





// computational routines for 2x2 matrices


template<typename T>
inline T det(const Mat2x2<T>& A)
{
    return A.v11 * A.v22 - A.v12 * A.v21;
}

template<typename T>
inline T det(const SMat2x2<T>& A)
{
    return A.v11 * A.v22 - A.v12 * A.v12;
}

template<typename T>
inline void inv(const Mat2x2<T>& A, Mat2x2<T>& R)
{
    T k = 1 / det(A);
    
    R.v11 = k * A.v22;
    R.v21 = - k * A.v21;
    R.v12 = - k * A.v12;
    R.v22 = k * A.v11;
}

template<typename T>
inline void inv(const SMat2x2<T>& A, SMat2x2<T>& R)
{
    T k = 1 / det(A);
    
    R.v11 = k * A.v22;
    R.v12 = -k * A.v12;
    R.v22 = k * A.v11;
}
        

template<typename T>
inline T trace(const Mat2x2<T>& A)
{
    return A.v11 + A.v22;
}

template<typename T>
inline T trace(const SMat2x2<T>& A)
{
    return A.v11 + A.v22;
}
        

template<typename T>
inline void polar_(T v11, T v12, T v22, T& a, T& b, T& theta)
{
    if (v12 == 0)
    {
        if (v11 >= v22)
        {
            a = v11;
            b = v22;
            theta = 0;
        }
        else
        {
            a = v22;
            b = v11;
            theta = M_PI_2;
        }
    }
    else
    {
        T y = 2 * v12;
        T x = v11 - v22;
        
        T d = std::sqrt(x * x + y * y);
        T s = v11 + v22;
        
        a = (s + d) / 2;
        b = (s - d) / 2;   
        
        theta = std::atan2(y, x) / 2;                
    }
}


template<typename T>
inline void polar(const SMat2x2<T>& A, Polar2<T>& R)
{
    polar_(A.v11, A.v12, A.v22, R.a, R.b, R.theta);        
}


template<typename T>
inline void polar(const Mat2x2<T>& A, Polar2<T>& R)
{
    polar_(A.v11, A.v12, A.v22, R.a, R.b, R.theta);        
}


template<typename T>
inline void polar2mat(const Polar2<T>& p, SMat2x2<T>& R)
{
    T c = std::cos(p.theta);
    T s = std::sin(p.theta);
    
    T c2 = c * c;
    T s2 = s * s;
    
    T a = p.a;
    T b = p.b;
    
    R.v11 = a * c2 + b * s2;
    R.v12 = (a - b) * c * s;
    R.v22 = a * s2 + b * c2;    
}


template<typename T>
inline void polar2mat(const Polar2<T>& p, Mat2x2<T>& R)
{
    T c = std::cos(p.theta);
    T s = std::sin(p.theta);
    
    T c2 = c * c;
    T s2 = s * s;
    
    T a = p.a;
    T b = p.b;
    
    R.v11 = a * c2 + b * s2;
    R.v12 = R.v21 = (a - b) * c * s;
    R.v22 = a * s2 + b * c2;    
}



template<typename T>
inline void sqrtm(const SMat2x2<T>& A, SMat2x2<T>& R)
{
    Polar2<T> p;
    polar(A, p);
    p.a = std::sqrt(p.a);
    p.b = std::sqrt(p.b);
    polar2mat(p, R);    
}

template<typename T>
inline void sqrtm(const Mat2x2<T>& A, Mat2x2<T>& R)
{
    Polar2<T> p;
    polar(A, p);
    p.a = std::sqrt(p.a);
    p.b = std::sqrt(p.b);
    polar2mat(p, R);    
}


template<typename T>
inline void chol_(T v11, T v21, T v22, T& r11, T& r21, T& r22)
{
    if (v11 == 0 && v21 == 0)
    {
        r11 = 0;
        r21 = 0;
        r22 = std::sqrt(v22);        
    }
    else
    {
        r11 = std::sqrt(v11);
        r21 = v21 / r11;
        r22 = std::sqrt(v22 - r21 * r21);
    }
}



template<typename T>
inline void chol(const Mat2x2<T>& A, Mat2x2<T>& R)
{
    chol_(A.v11, A.v21, A.v22, R.v11, R.v21, R.v22);
    R.v12 = 0;
}

template<typename T>
inline void chol(const SMat2x2<T>& A, LTMat2x2<T>& R)
{
    chol_(A.v11, A.v12, A.v22, R.v11, R.v21, R.v22);
}





