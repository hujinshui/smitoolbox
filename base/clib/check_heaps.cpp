// A program for checking the correctness of heap implementation


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "array.h"
#include "heaps.h"

using namespace smi;

enum Action
{
    ADD_KEY = 0,
    REMOVE_ROOT = 1,
    SET_KEY = 2    
};


Action pick_action()
{
    int a = rand() % 3;
    return (Action)a;
}

double pick_key()
{
    return (double)rand() / RAND_MAX;    
}

size_t pick_index(size_t maxIdx)
{
    return (size_t)rand() % (maxIdx + 1);
}



template<typename T, typename TOrd>
void verify_heap(const BinaryHeap<T, TOrd>& H)
{
    typedef typename BinaryHeap<T, TOrd>::btree_type btree_t;
    typedef typename BinaryHeap<T, TOrd>::trnode node_t;
    
    const btree_t& bt = H.btree();
    TOrd ord;
    
    // verify size
    
    if (H.size() != bt.size())
    {
        printf("heap size != tree size\n");
        exit(1);
    }
        
    // verify heap order condition
    
    if (!bt.empty())
    {
        for (int ip = 1; ip <= bt.size(); ++ip)
        {
            node_t p = ip;
            node_t par = p.parent();
            
            if (!H.is_inheap(bt[p]))
            {
                printf("a node in btree is tagged as not in heap\n");
                exit(1);
            }                        
            
            if (!par.is_null())
            {
                if (ord( H.key_at_node(p), H.key_at_node(par) ))
                {
                    printf("heap condition violated \n");
                    exit(1);
                }
            }  
            
            size_t idx = bt[p];
            if (p != H.get_node_by_index(idx))
            {
                printf("Back link (idx -> node) inconsistency\n");
                exit(1);
            }
        }
    } 
    
    printf("verified heap [size = %d]\n", H.size());
}




template<class THeap>
class TestContext
{
public:
    TestContext() : m_max_index(-1) { }    
    
    void random_action()
    {
        do_action(pick_action());
    }    
    
    void do_action(Action a)
    {
        switch (a)
        {
            case ADD_KEY:
                do_add_key();
                break;
            case REMOVE_ROOT:
                do_remove_root();
                break;
            case SET_KEY:
                do_set_key();
                break;
        }
    }
    
    
    void do_make_heap(size_t n0)
    {
        printf("To make heap with n0 = %d\n", n0);
        
        Array<double> src(n0);
        for (size_t i = 0; i < n0; ++i) src[i] = pick_key();
        
        m_heap.make_heap(src.begin(), src.end());
        _verify_heap();
        
        m_max_index = n0-1;
    }
    
    void do_add_key()
    {
        double kv = pick_key();
        
        printf("To add key %.4f ...\n", kv);
        size_t new_idx = m_heap.add_key(kv);
        
        if (new_idx == (m_max_index + 1))
        {
            m_max_index = new_idx;
        }
        else
        {
            printf("max index rule violated.\n");
            exit(1);
        }  
        
        _verify_heap();
    }
    
    void do_remove_root()
    {        
        if (m_heap.empty()) return;
        
        printf("To delete root ...\n");         
        m_heap.delete_root();
        
        _verify_heap();        
    }
    
    void do_set_key()
    {        
        int i = pick_index(m_max_index);
        
        if (m_heap.is_inheap(i))
        {
            double kv = pick_key();
            printf("To set key [%d] = %.4f ...\n", i, kv);
            m_heap.set_key(i, kv);
            
            _verify_heap();
        }              
    }    
    
    
private:
    void _verify_heap()
    {
        typedef typename THeap::key_type key_t;
        typedef typename THeap::key_order ord_t;
        
        verify_heap<key_t, ord_t>(m_heap);
    }    
    
private:
    THeap m_heap;
    size_t m_max_index;
};




int main(int argc, char *argv[])
{
    srand ( time(NULL) );
    
    typedef BinaryHeap<double> BHeap;
    
    if (argc != 3)
    {
        printf("Require exactly two input arguments.\n");
        exit(1);
    }
    
    size_t n0 = (size_t)atoi(argv[1]);
    size_t na = (size_t)atoi(argv[2]);
        
    if (n0 == 0 || na == 0)
    {
        printf("Both inputs should be positive integers.\n");
        exit(1);
    }
        
    printf("Test on n0 = %d with %d actions ......\n", n0, na);
    
    TestContext<BHeap> ctx;
    
    ctx.do_make_heap(n0);
         
    for (size_t i = 0; i < na; ++i)
    {
        ctx.random_action();
    }
    
    printf("\n");
}


