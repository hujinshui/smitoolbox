/********************************************************************
 *
 *  C++ mex implementation of pairwise hamming distance computing
 *
 *  Created by Dahua Lin, on Jun 3, 2008
 *
 ********************************************************************/

#include <bcslib/matlab/bcs_mex.h>

using namespace bcs;
using namespace bcs::matlab;


template<typename T>
double hamdist(const const_aview1d<T>& x, const const_aview1d<T>& y)
{
    index_t m = (index_t)x.nelems();
    
    unsigned int s = 0;
    for (index_t i = 0; i < m; ++i)
    {
        if (x[i] != y[i]) ++ s; 
    }
    
    return s;
}
         
        

template<typename T>
marray do_pwhamdist(const_marray mX1, const_marray mX2)
{
    size_t m = mX1.nrows();
    size_t n1 = mX1.ncolumns();
    size_t n2 = mX2.ncolumns();
    
    const_aview2d<T, column_major_t> X1 = view2d<T>(mX1);
    const_aview2d<T, column_major_t> X2 = view2d<T>(mX2);
    
    marray mD = create_marray<double>(n1, n2);
    aview2d<double, column_major_t> D = view2d<double>(mD);
    
    for (size_t j = 0; j < n2; ++j)
    {
        const_aview1d<T> vj = X2.column(j);
        
        for (size_t i = 0; i < n1; ++i)
        {
            const_aview1d<T> vi = X1.column(i);
            D(i, j) = hamdist(vi, vj);
        }
    }
    
    return mD;
}


/*
 * entry
 *
 * Input
 *    - X1, X2: the input matrices d x n1, and d x n2
 *
 * Output
 *    - dists:  n1 x n2 pairwise distance matrix
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    // take inputs
    
    const_marray mX1(prhs[0]);
    const_marray mX2(prhs[1]);
    
    if (!( mX1.is_matrix() && (mX1.is_numeric() || mX1.is_logical()) && !mX1.is_sparse() ))  
    {
        throw mexception("pwhamdist:invalidarg", 
            "X1 should be a non-sparse matrix.");    
    }
    
    if (!( mX2.is_matrix() && (mX2.is_numeric() || mX2.is_logical()) && !mX2.is_sparse() ))    
    {
        throw mexception("pwhamdist:invalidarg", 
            "X2 should be a non-sparse matrix.");
    }
    
    if ( mX1.nrows() != mX2.nrows() )
    {
        throw mexception("pwhamdist:invalidarg", 
            "X1 and X2 should have the same number of rows.");
    }
    
    if ( mX1.class_id() != mX2.class_id() )
    {
        throw mexception("pwhamdist:invalidarg", 
            "X1 and X2 should have the same value type.");
    }
    
    // main delegation
    
    marray mD;
    
    switch (mX1.class_id())
    {
        case mxDOUBLE_CLASS:
            mD = do_pwhamdist<double>(mX1, mX2);
            break;
            
        case mxSINGLE_CLASS:
            mD = do_pwhamdist<float>(mX1, mX2);
            break;
            
        case mxLOGICAL_CLASS:
            mD = do_pwhamdist<bool>(mX1, mX2);
            break;
            
        case mxINT32_CLASS:
            mD = do_pwhamdist<int32_t>(mX1, mX2);
            break;
            
        case mxUINT32_CLASS:
            mD = do_pwhamdist<uint32_t>(mX1, mX2);
            break;       
            
        case mxINT16_CLASS:
            mD = do_pwhamdist<int16_t>(mX1, mX2);
            break;
            
        case mxUINT16_CLASS:
            mD = do_pwhamdist<uint16_t>(mX1, mX2);
            break;
            
        case mxINT8_CLASS:
            mD = do_pwhamdist<int8_t>(mX1, mX2);
            break;
            
        case mxUINT8_CLASS:
            mD = do_pwhamdist<uint8_t>(mX1, mX2);
            break;             
    }
    
    
    // output
    
    plhs[0] = mD.mx_ptr();
}


BCSMEX_MAINDEF

