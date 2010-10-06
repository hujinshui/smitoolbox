/********************************************************************
 *
 *  clib_ex2.cpp
 *
 *  An example showing how to use heaps in clib
 *
 ********************************************************************/

#define MATLAB_INSPECT_HEAPS
#define MATLAB_DUMP_HEAP_ACTION

#include "marray.h"
#include "heaps.h"

using namespace smi;

typedef double key_type;


template class CompleteBinaryTree<int>;
template class BinaryHeap<key_type>;

typedef BinaryHeap<key_type> BHeap;


inline void print_node(const ConstArray<key_type>& src, int i)
{
    mexPrintf("%g", double(src[i]));  
}


void print_array(const char *title, const ConstArray<key_type>& src)
{
    mexPrintf(title);
    for (int i = 0; i < src.size(); ++i)
    {
        print_node(src, i);
        mexPrintf(" ");
    }
    mexPrintf("\n");
}


template<class THeap>
void print_heap(const THeap& h)
{
    mexPrintf("heap:\n");
    h.dump_heap_to_matlab();
    mexPrintf("imap:\n");
    h.dump_imap_to_matlab();
}


template<class THeap>
void run_on(const ConstArray<key_type>& src)
{
    int n = src.size();
    Array<key_type> buffer(n);
    
    print_array("source:\n", src);
    mexPrintf("\n");
    
    THeap heap1(2 * n);
        
    mexPrintf("Make heap:\n");
    heap1.make_heap(n, src.data());
    print_heap(heap1);
    mexPrintf("\n");
    
    mexPrintf("Extract all:\n");
    heap_extract_n(heap1, heap1.size(), buffer.data());
    print_heap(heap1);
    print_array("result:\n", buffer);
    mexPrintf("\n");
    
    THeap& heap2 = heap1;
    
    mexPrintf("Refill heap by sequential insertion:\n");
    for (int i = 0; i < n; ++i)
    {
        mexPrintf("add ");
        print_node(src, i);
        mexPrintf("\n");
        heap2.add_key(src[i]);
    }
    print_heap(heap2);
    mexPrintf("\n");    
}


/**
 * Main entry:
 *
 * Input:
 *   [0]: src  the source array (should be a double vector)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgTxt("The number of inputs should be exactly 1.");
    
    MArray mSrc(prhs[0]);
    
    if (!(mSrc.is_double() && !mSrc.is_sparse()))
    {
        mexErrMsgTxt("The src should be a non-sparse double vector.");
    }
    
    ConstArray<key_type> src(mSrc.nelems(), mSrc.get_data<key_type>());
    
    run_on<BHeap>(src);
}

