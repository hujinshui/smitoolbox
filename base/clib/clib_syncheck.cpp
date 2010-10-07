// a special file for syntax checking of template libraries


#include "matlab_types.h"
#include "marray.h"

#include "array.h"
#include "bintree.h"
#include "heaps.h"

using namespace smi;

template class Array<double>;
template class Array<bool>;

template class CompleteBinaryTree<double>;
template class BinaryHeap<double>;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
}