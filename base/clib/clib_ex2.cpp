/********************************************************************
 *
 *  clib_ex2.cpp
 *
 *  An example showing how to use heaps in clib
 *
 ********************************************************************/

#define SMI_BINARY_HEAP_INSPECTION

#include "marray.h"
#include "array.h"
#include "heaps.h"

using namespace smi;


typedef BinaryHeap<double> BHeap;

template<typename TArray>
void print_array(const char* title, const TArray& a)
{
    mexPrintf("%s ", title);
    for (typename TArray::const_iterator it = a.begin(); it != a.end(); ++it)
    {
        mexPrintf("%g ", double(*it));
    }
    mexPrintf("\n");
}


class Monitor : public BHeap::IMonitor
{
public:
    typedef BHeap::index_type index_type;
    typedef BHeap::key_type key_type;
    
public:
    virtual void on_node_swaped(const BHeap& H, index_type i, index_type j)
    {
        mexPrintf("swap %g <--> %g\n", H.get_key(i), H.get_key(j));
    }
        
    virtual void on_node_appended(const BHeap& H, index_type i)
    {
        mexPrintf("append [%d] = %g\n", i, H.get_key(i));
    }

    virtual void on_upheaping(const BHeap& H, index_type i)
    {
        mexPrintf("To upheap %g\n", H.get_key(i));
    }

    virtual void on_downheaping(const BHeap& H, index_type i)
    {
        mexPrintf("To downheap %g\n", H.get_key(i));
    }

    virtual void on_key_setting(const BHeap& H, index_type i, const key_type& kv)
    {        
        mexPrintf("To set [%d] = %g\n", i, kv);
    }

    virtual void on_key_adding(const BHeap& H, const key_type& kv)
    {
        mexPrintf("to add %g\n", kv);
    }

    virtual void on_heap_making(const BHeap& H)
    {
        mexPrintf("to make-heap\n");
    }

    virtual void on_heap_made(const BHeap& H)
    {
        mexPrintf("heap-made => ");
        print_heap(H);
    }

    virtual void on_clearing(const BHeap& H)
    {
        mexPrintf("to clear-heap\n");
    }

    virtual void on_cleared(const BHeap& H)
    {
        mexPrintf("heap-cleared => ");
        print_heap(H);
    }

    virtual void on_root_deleting(const BHeap& H)
    {
        mexPrintf("to delete root %g\n", H.root_key());
    }

    virtual void on_root_deleted(const BHeap& H)
    {
        mexPrintf("root-deleted => ");
        print_heap(H);
    }
    
    void print_heap(const BHeap& H)
    {
        const BHeap::btree_type& T = H.btree(); 
        
        for (BHeap::btree_type::const_iterator it = T.begin(); it != T.end(); ++it)
        {
            mexPrintf("%g ", H.get_key(*it));
        }
        mexPrintf("\n");
    }
};





template<class THeap>
void run_on(const CRefArray<double>& src)
{
    std::vector<double> results;
    
    // print_array("Source:", src);
    
    Monitor mon;
    
    THeap heap1;        
    heap1.set_monitor(&mon);
        
    mexPrintf("Test 1: Make Heap\n");
    mexPrintf("===========================\n");
    heap1.make_heap(src.begin(), src.end());
    mexPrintf("\n");
    
    mexPrintf("Test 2: Heap sort\n");
    mexPrintf("===========================\n");
    while (!heap1.empty())
    {
        results.push_back(heap1.root_key());
        heap1.delete_root();
    }
    mexPrintf("----------------------\n");
    print_array("result:", results);
    mexPrintf("\n");
    
    THeap& heap2 = heap1;
    
    mexPrintf("Test 3: Re-insertion\n");
    mexPrintf("===========================\n");
    for (size_t i = 0; i < src.size(); ++i)
    {
        heap2.add_key(src[i]);
    }
    mexPrintf("----------------------\n");
    mexPrintf("heap after insertion: ");
    mon.print_heap(heap2);    
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
    
    CRefArray<double> src(mSrc.nelems(), mSrc.get_data<double>());
    
    run_on<BHeap>(src);
}

