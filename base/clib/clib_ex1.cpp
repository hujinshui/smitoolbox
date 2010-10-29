/********************************************************************
 *
 *  clib_ex1.cpp
 *
 *  An example showing how to use clib
 *
 ********************************************************************/

#include <mex.h>

#include "vector.h"
#include "matrix.h"
#include "marray.h"

using namespace smi;


template class VectorCView<double>;
template class VectorView<double>;
template class MultiVectorCView<double>;
template class MultiVectorView<double>;
template class MatrixCView<double>;
template class MatrixView<double>;

void print_matrix(const MatrixCView<double>& M)
{
    for (int i = 0; i < M.nrows(); ++i)
    {
        for (int j = 0; j < M.ncols(); ++j)
        {
            mexPrintf("%.4g ", M(i, j));
        }
        mexPrintf("\n");
    }
}

void print_vector(const VectorCView<double>& v)
{
    for (int i = 0; i < v.nelems(); ++i)
    {
        mexPrintf("%.4g ", v(i));
    }    
}

void print_vectors(const MultiVectorCView<double>& V)
{        
    for (int i = 0; i < V.nvecs(); ++i)
    {
        mexPrintf("[%d]: ", i);
        print_vector(V[i]);
        mexPrintf("\n");
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgTxt("The number of inputs should be 1.");
    
    MArray mA(prhs[0]);
    
    if (!(mA.is_real_nonsparse_matrix() && mA.is_double()))
    {
        mexErrMsgTxt("The input array should be a real double matrix.");
    }
    
    MatrixCView<double> A = mA.to_matrix<double>();
    
    mexPrintf("Matrix:\n");
    print_matrix(A);
    mexPrintf("\n");
    
    mexPrintf("Flat:\n");
    print_vector(A.flat());
    mexPrintf("\n\n");
    
    mexPrintf("Rows:\n");
    print_vectors(A.rows());
    mexPrintf("\n");
    
    mexPrintf("Columns:\n");
    print_vectors(A.columns());        
    mexPrintf("\n");
}




