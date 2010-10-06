/********************************************************************
 *
 *  sorting.h
 *
 *  All sorts of sorting algorithm
 *
 *  Created by Dahua Lin, on Oct 5, 2010
 *
 ********************************************************************/


#ifndef SMI_CLIB_SORTING_H
#define SMI_CLIB_SORTING_H


// To dump steps of sorts
//
// define a macro  MATLAB_DUMP_QSORT_ACTION
//

#if defined(MATLAB_DUMP_QSORT_ACTION)
#include <mex.h>    
#endif

#include "data_structs.h"
#include "heaps.h"

#include <math.h>


namespace smi
{
        
    
/************************************************
 *
 *  Auxiliary structures
 *
 ************************************************/    
    
template<typename T, bool OutputIndex>    
struct SortEntry;    
    
    
template<typename T>
struct SortEntry<T, false>
{
    // fields
    
    T *data;
    int n;
    
    // constructors
    
    SortEntry() { }    
    SortEntry(T *data_, int n_) : data(data_), n(n_) { }
    
    // operations
    
    int left_index() const { return 0; }    
    int right_index() const { return n-1; }    
    int mid_index() const { return n >> 1; }
    
    void swap(int i, int j)
    {
        T t = data[i];
        data[i] = data[j];
        data[j] = t;
    }    
    
    SortEntry<T, false> sub_entry(int ifirst, int ilast) const
    {
        return SortEntry<T, false>(data + ifirst, ilast - ifirst + 1);
    }
};


template<typename T>
struct SortEntry<T, true>
{
    // fields
    
    T *data;
    int *indices;
    int n;
    
    // constructors
    
    SortEntry() { }    
    SortEntry(T *data_, int *inds, int n_) : data(data_), indices(inds), n(n_) { }
    
    // operations
    
    int left_index() const { return 0; }    
    int right_index() const { return n-1; }
    int mid_index() const { return n >> 1; }

    
    void swap(int i, int j)
    {
        T t = data[i];
        data[i] = data[j];
        data[j] = t;
        
        int ti = indices[i];
        indices[i] = indices[j];
        indices[j] = ti;
    }
    
    SortEntry<T, true> sub_entry(int ifirst, int ilast) const
    {
        return SortEntry<T, true>(data + ifirst, indices + ifirst, ilast - ifirst + 1);
    }
};


template<typename T>
inline SortEntry<T, false> make_sort_entry(T *data, int n)
{
    return SortEntry<T, false>(data, n);
}

template<typename T>
inline SortEntry<T, true> make_sort_entry(T *data, int *indices, int n)
{
    return SortEntry<T, true>(data, indices, n);
}


template<typename T, typename TOrd>
inline bool is_sorted(const T *data, int n, TOrd ord)
{
    if (n > 1)
    {
        for (int i = 1; i < n; ++i)
        {
            if (ord(data[i], data[i-1])) return false;
        }
        return true;
    }
    else
    {
        return true;
    }        
}

template<typename T>
inline bool is_sorted_asc(const T *data, int n)
{
    return is_sorted(data, n, less<T>());
}

template<typename T>
inline bool is_sorted_desc(const T *data, int n)
{
    return is_sorted(data, n, greater<T>());
}




    
/************************************************
 *
 *  Heap sort
 *
 ************************************************/       
    
      
template<typename T, typename TOrd>    
void heapsort(int n, const T *src, T* val_out, int *idx_out)
{
    BinaryHeap<T, TOrd> heap(n);
    
    heap.make_heap(n, src);
    
    for (int t = 0; t < n; ++t)
    {
        int i;
        *(val_out++) = heap.extract_root(i);
        *(idx_out++) = i;
    }
}


template<typename T>
inline void heapsort_asc(int n, const T *src, T *val_out, int *idx_out)        
{
    heapsort<T, less<T> >(n, src, val_out, idx_out);
}


template<typename T>
inline void heapsort_desc(int n, const T *src, T *val_out, int *idx_out)        
{
    heapsort<T, greater<T> >(n, src, val_out, idx_out);
}



/************************************************
 *
 *  Quick sort
 *
 ************************************************/


template<typename T, bool OutputIndex, typename TOrd>
inline void sort2(SortEntry<T, OutputIndex> e, TOrd ord)
{
    if (ord(e.data[1], e.data[0]))
    {
        e.swap(0, 1);
    }
}

template<typename T, bool OutputIndex, typename TOrd>
inline void sort3(SortEntry<T, OutputIndex> e, TOrd ord)
{
    if (ord(e.data[1], e.data[0])) e.swap(0, 1);        
    if (ord(e.data[2], e.data[1])) e.swap(1, 2);        
    if (ord(e.data[1], e.data[0])) e.swap(0, 1);    
}


template<typename T, bool OutputIndex, typename TOrd>
class QuickSortPartition
{
public:
    
    typedef SortEntry<T, OutputIndex> entry_type;
    
    QuickSortPartition(entry_type e, TOrd ord_): m_e(e), ord(ord_), m_done(false)    
    {        
    }
    
    bool done() const
    {
        return m_done;
    }
                    
    
    entry_type left_entry() const  // only use this when further recursion is needed
    {
        return m_e.sub_entry(m_e.left_index(), m_pi-1);        
    }
    
    entry_type right_entry() const  // only use this when further recursion is needed
    {
        return m_e.sub_entry(m_pi+1, m_e.right_index());
    }
    
    
    /**
     * Performs partitioning
     * 
     * @return whether further recursion is needed
     */
    bool run()
    {
        #ifdef MATLAB_DUMP_QSORT_ACTION
                mdump_progress("On partition: ", -1, -1);
        #endif  
        
        initialize();
        
        if (!done())
        {                                  
            T pv = m_pv;
            
            // swap pivot to last
            
            int last_i = m_e.right_index();                        
            m_e.swap(m_pi, last_i);
            
            #ifdef MATLAB_DUMP_QSORT_ACTION
                    mdump_progress("swap pivot: ", 0, 0);
            #endif
            
            // scan each element
            
            int si = 0;            
            for (int i = 0; i < last_i; ++i)
            {
                T v = m_e.data[i];
                
                if (!ord(pv, v))  // v --> left
                {
                    if (i > si)
                    {
                        m_e.swap(si, i);
                    }
                    ++ si;
                }
                
                #ifdef MATLAB_DUMP_QSORT_ACTION
                    mdump_progress("scan next: ", si, i+1);
                #endif
            }
            
            // swap back the pivot
            
            m_pi = si;
            m_e.swap(si++, last_i);              
            
            #ifdef MATLAB_DUMP_QSORT_ACTION
                    mdump_progress("pivot back: ", si, -1);
            #endif
            
            m_done = true;
            
            #ifdef MATLAB_DUMP_QSORT_ACTION
                    mexPrintf("\n");
            #endif
                    
            return true;
        }
        else
        {
            #ifdef MATLAB_DUMP_QSORT_ACTION
                    mexPrintf("\n");
            #endif
                    
            return false;
        }
    }
    
    
private:
        
    void initialize()
    {                
        if (m_e.n == 2)
        {
            sort2(m_e, ord);                          
            m_done = true;
            
            #ifdef MATLAB_DUMP_QSORT_ACTION
                    mdump_progress("sort2: ", -1, -1);
            #endif
        }
        else if (m_e.n == 3)
        {
            sort3(m_e, ord);            
            m_done = true;
            
            #ifdef MATLAB_DUMP_QSORT_ACTION
                    mdump_progress("sort 3: ", -1, -1);
            #endif
        }
        else if (m_e.n > 3)
        {
            choose_pivot();            
            m_done = false;
        }
        else
        {
            m_done = true;
        }
    }
    
    void choose_pivot()
    {
        T v1 = m_e.data[m_e.left_index()];
        T v2 = m_e.data[m_e.mid_index()];
        T v3 = m_e.data[m_e.right_index()];
        
        if (ord(v1, v3))
        {
            if (ord(v2, v3))
            {
                if (ord(v1, v2))  // v1 < v2 < v3
                { 
                    m_pi = m_e.mid_index(); m_pv = v2; 
                }
                else  // v2 <= v1 < v3
                { 
                    m_pi = m_e.left_index(); m_pv = v1; 
                }
            }
            else // v1 < v3 <= v2
            {
                m_pi = m_e.right_index(); m_pv = v3;
            }
        }
        else // v1 >= v3
        {
            if (ord(v2, v1))
            {
                if (ord(v3, v2)) // v3 < v2 < v1
                { 
                    m_pi = m_e.mid_index(); m_pv = v2; 
                }
                else  // v2 <= v3 <= v1
                {
                    m_pi = m_e.right_index(); m_pv = v3;
                }
            }
            else // v3 <= v1 <= v2
            {
                m_pi = m_e.left_index(); m_pv = v1;
            }
        }
    }    
       
#ifdef MATLAB_DUMP_QSORT_ACTION

private:
    void mdump_progress(const char *title, int si, int fi)
    {
        mexPrintf(title);
        
        for (int i = 0; i < m_e.right_index(); ++i)
        {
            if (i == si)
            {
                mexPrintf("| ");                                
            }
            
            if (i == fi)
            {
                mexPrintf("* ");
            }
            
            mexPrintf("%.4g ", double(m_e.data[i]));                        
        }
        
        if (fi >= 0)
        {
            mexPrintf("(%.4g) ", double(m_e.data[m_e.right_index()]));
        }
        else
        {
            mexPrintf("%.4g ", double(m_e.data[m_e.right_index()]));
        }
        
        mexPrintf("\n");
    }
    
#endif
    
    
public:
    entry_type m_e;
    TOrd ord;    
    bool m_done;
    
    int m_pi;
    T m_pv;
        
}; // end QuickSortIterator



template<typename T, bool OutputIndex, typename TOrd>
void quicksort_ip(SortEntry<T, OutputIndex> e, TOrd ord)
{
    typedef SortEntry<T, OutputIndex> entry_t;
    typedef QuickSortPartition<T, OutputIndex, TOrd> par_t;
    
    
    int max_stack_size = (int)ceil(log2((double)e.n));     
    Stack<entry_t> stack(e.n+1);
    
    stack.push(e);
    mexPrintf("push %d => %d\n", e.n, stack.size());
    
    int ass = 0;
    
    while (!stack.empty())
    {
        entry_t ce = stack.top();
        stack.pop();          
        mexPrintf("pop %d => %d\n", ce.n, stack.size());
                
        bool b = true; // whether to continue inner loop
        while (b)
        {
            b = false;
            par_t par(ce, ord);
            
            if (par.run())
            {
                entry_t el = par.left_entry();
                entry_t er = par.right_entry();
            
                if (el.n >= er.n)
                {
                    if (el.n > 1) 
                    {
                        ce = el;
                        b = true;
                        
                        if (er.n > 1) stack.push(er);
                        mexPrintf("push %d => %d\n", er.n, stack.size());
                    }                                        
                }
                else
                {
                    if (er.n > 1)
                    {
                        ce = er;
                        b = true;
                        
                        if (el.n > 1) stack.push(el);
                        mexPrintf("push %d => %d\n", el.n, stack.size());
                    }                    
                }
            } 
            
            if (stack.size() > ass) ass = stack.size();
            
        } // inner loop        
    }
    
    mexPrintf("ass = %d\n", ass);
    mexPrintf("max_stack_size = %d\n", max_stack_size);
    
}



template<typename T, typename TOrd>
void quicksort(int n, const T *src, T* val_out, int* idx_out = 0)
{   
    if (n >= 1)
    {
        // do copying if necessary
    
        if (val_out > src)
        {
            for (int i = n-1; i >= 0; --i) val_out[i] = src[i];
        }
        else if (val_out < src)
        {
            for (int i = 0; i < n; ++i) val_out[i] = src[i];
        }

        // main 

        if (idx_out == 0)
        {   
            if (n > 1)
                quicksort_ip(make_sort_entry(val_out, n), TOrd());
        }
        else
        {
            if (n > 1)
            {
                for (int i = 0; i < n; ++i) idx_out[i] = i;
                quicksort_ip(make_sort_entry(val_out, idx_out, n), TOrd());
            }
            else
            {
                *idx_out = 0;
            }
        }
    }
}


template<typename T>
inline void quicksort_asc(int n, const T *src, T* val_out, int* idx_out = 0)
{
    quicksort<T, less<T> >(n, src, val_out, idx_out);
}

template<typename T>
inline void quicksort_desc(int n, const T *src, T* val_out, int* idx_out = 0)
{
    quicksort<T, greater<T> >(n, src, val_out, idx_out);
}


    
}


#endif





