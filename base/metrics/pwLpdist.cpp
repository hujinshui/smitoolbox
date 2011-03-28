/********************************************************************
 *
 *  The C++ mex implementation of pairwise Minkowski distance,
 *  including L1-distance, L2-distance, and L-inf distance, as
 *  well as generic Minkowski distance.
 *
 *  Created by Dahua Lin, on Jun 2, 2008
 *
 ********************************************************************/


#include <bcslib/matlab/bcs_mex.h>
#include <bcslib/veccomp/veccalc.h>
#include <bcslib/array/array_eval.h>

using namespace bcs;
using namespace bcs::matlab;


template<typename T>
T get_power(const_marray mP, bool &is_finite)
{
    if ( !(mP.is_numeric() && mP.is_scalar() && !mP.is_sparse()) )
    {
        throw mexception("pwLpdist:invalidarg", 
            "p should be a numeric scalar.");
    }
    
    T p;    
    switch (mP.class_id())
    {
        case mxDOUBLE_CLASS:
            p = (T)mP.get_scalar<double>();
            break;
        case mxSINGLE_CLASS:
            p = (T)mP.get_scalar<float>();
            break;
        case mxINT32_CLASS:
            p = (T)mP.get_scalar<int32_t>();
            break;
        case mxUINT32_CLASS:
            p = (T)mP.get_scalar<uint32_t>();
            break;
        default:
            throw mexception("pwLpdist:invalidarg", 
                "p should be of type double, single, int32, or uint32.");                    
    }
    
    if (p < 1)
    {
        throw mexception("pwLpdist:invalidarg", 
            "p should be a positive value no less than 1.");
    }
    
    is_finite = mxIsFinite(p);    
}


template<typename T>
marray do_pwLpdist(const_marray mX1, const_marray mX2, const T& p, bool is_finite)
{
    const_aview2d<T, column_major_t> X1 = view2d<T>(mX1);
    const_aview2d<T, column_major_t> X2 = view2d<T>(mX2);
    
    size_t m = X1.nrows();
    size_t n1 = X1.ncolumns();
    size_t n2 = X2.ncolumns();
    
    array1d<T> temp(m);
    
    marray mD = create_marray<T>(n1, n2);
    aview2d<T, column_major_t> D = view2d<T>(mD);
        
    for (index_t j = 0; j < (index_t)n2; ++j)
    {
        const_aview1d<T> vj = X2.column(j);
        
        for (index_t i = 0; i < (index_t)n1; ++i)
        {
            const_aview1d<T> vi = X1.column(i);
            
            vec_sub(m, vi.pbase(), vj.pbase(), temp.pbase());
            
            if (is_finite)
            {
                D(i, j) = Lpnorm(temp, p);
            }
            else
            {
                D(i, j) = Linfnorm(temp);
            }            
        }
    }        
    
    return mD;
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
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    // take inputs
    
    const_marray mX1(prhs[0]);
    const_marray mX2(prhs[1]);
    const_marray mP(prhs[2]);
    
    if (!( mX1.is_matrix() && mX1.is_float() && !mX1.is_sparse() ))  
    {
        throw mexception("pwLpdist:invalidarg", 
            "X1 should be a non-sparse float-typed matrix.");    
    }
    
    if (!( mX2.is_matrix() && mX2.is_float() && !mX2.is_sparse() ))    
    {
        throw mexception("pwLpdist:invalidarg", 
            "X2 should be a float-typed matrix.");
    }
    
    if ( !( mX1.nrows() == mX2.nrows() ))
    {
        throw mexception("pwLpdist:invalidarg", 
            "X1 and X2 should have the same number of rows.");
    }
    
    // main delegation
    
    marray mD;
    bool is_finite;
    
    if (mX1.is_double() && mX2.is_double())
    {        
        double p = get_power<double>(mP, is_finite);
        mD = do_pwLpdist<double>(mX1, mX2, p, is_finite);
    }
    else if (mX1.is_single() && mX2.is_single())
    {
        float p = get_power<float>(mP, is_finite);
        mD = do_pwLpdist<double>(mX1, mX2, p, is_finite);
    }
    else
    {
        throw mexception("pwLpdist:invalidarg", 
            "X1 and X2 should have the same value type.");
    }
    
    // output
    
    plhs[0] = mD.mx_ptr();
}


BCSMEX_MAINDEF
