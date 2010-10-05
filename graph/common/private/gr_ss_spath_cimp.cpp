/********************************************************************
 *
 *  gr_ss_spath_cimp.cpp
 *
 *  The C++ mex implementation for several single-source 
 *  shortest path solver
 *
 *  Created by Dahua Lin, on Oct 4, 2010
 *
 ********************************************************************/

#include "../../clib/mgraph.h"
#include "../../clib/graph_spath.h"

using namespace smi;

template<class TAlgIter, typename T>
struct base_getter
{
    typedef T value_type;
    
    const TAlgIter& m_it;
    const SeqList<int>& m_seq;
    
    base_getter(const TAlgIter& it) : m_it(it), m_seq(it.close_sequence()) { }
    
    int count() const { return m_it.num_closed(); }
};

template<class TAlgIter>
struct vnode_getter : public base_getter<TAlgIter, int>
{
    vnode_getter(const TAlgIter& it) : base_getter<TAlgIter, int>(it) { }
    
    int operator() (int i) const { return this->m_seq[i] + 1; }
};

template<class TAlgIter, class TWeight>
struct sdist_getter : public base_getter<TAlgIter, TWeight>
{
    sdist_getter(const TAlgIter& it) : base_getter<TAlgIter, TWeight>(it) { }
    
    TWeight operator() (int i) const { return this->m_it.distance_of(this->m_seq[i]); }
};

template<class TAlgIter>
struct vpred_getter : public base_getter<TAlgIter, int>
{
    vpred_getter(const TAlgIter& it) : base_getter<TAlgIter, int>(it) { }
    
    int operator() (int i) const { return this->m_it.predecessor_of(this->m_seq[i]) + 1; }
};



template<typename TAlgIter, typename TWeight>
void extract_output(const TAlgIter& alg_it, 
        mxArray*& mxVs, mxArray*& mxDists, mxArray*& mxPreds)
{
    mxVs = iter_to_matlab_row(vnode_getter<TAlgIter>(alg_it));
    mxDists = iter_to_matlab_row(sdist_getter<TAlgIter, TWeight>(alg_it));
    mxPreds = iter_to_matlab_row(vpred_getter<TAlgIter>(alg_it));
}


template<typename T>
void do_solve(const MArray& mG, int s, char code, 
       mxArray*& mxVs, mxArray*& mxDists, mxArray*& mxPreds)
{
    RefWGraph<T> G = to_refwgraph<T>(mG);
    WAdjList<T> adjList = WAdjList<T>(G);

    if (code == 'a')
    {
        DAG_SPathIterator<T> sp_it(adjList);
        if (sp_it.initialize(s))
        {            
            sp_it.solve_all();        
            extract_output<DAG_SPathIterator<T>, T>(sp_it, mxVs, mxDists, mxPreds);
        }
        else
        {
            mexErrMsgIdAndTxt("gr_dag_spath:rterror", 
                    "The input graph is not an acyclic graph (the part accessible from s).");
        } 
    }
       
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    MArray mG(prhs[0]);
    MArray mS(prhs[1]);
    MArray mCode(prhs[2]);
               
    int s = mS.get_scalar<int>() - 1;
    char code = (char)(mCode.get_scalar<mxChar>());
  
    // main
    
    switch (get_graph_weight_class(mG))
    {
        case mxDOUBLE_CLASS:
            do_solve<double>(mG, s, code, plhs[0], plhs[1], plhs[2]);
            break;
        case mxSINGLE_CLASS:
            do_solve<float>(mG, s, code, plhs[0], plhs[1], plhs[2]);
            break;
        case mxINT32_CLASS:
            do_solve<int>(mG, s, code, plhs[0], plhs[1], plhs[2]);
            break;
        default:
            mexErrMsgIdAndTxt("gr_dag_spath:invalidarg", 
                    "The type of the weights must be either double, single, or int32.");
    }  
}

