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
    const TAlgIter& m_it;    
    base_getter(const TAlgIter& it) : m_it(it) { }
};

template<class TAlgIter>
struct vnode_getter : public base_getter<TAlgIter, int>
{
    vnode_getter(const TAlgIter& it) : base_getter<TAlgIter, int>(it) { }
    
    int operator() (int v) const 
    { 
        return v + 1; 
    }
};

template<class TAlgIter, class TWeight>
struct sdist_getter : public base_getter<TAlgIter, TWeight>
{
    sdist_getter(const TAlgIter& it) : base_getter<TAlgIter, TWeight>(it) { }    
    TWeight operator() (int v) const 
    { 
        return this->m_it.distance_of(v); 
    }
};

template<class TAlgIter>
struct vpred_getter : public base_getter<TAlgIter, int>
{
    vpred_getter(const TAlgIter& it) : base_getter<TAlgIter, int>(it) { }
    
    int operator() (int v) const 
    { 
        return this->m_it.predecessor_of(v) + 1; 
    }
};



template<typename TAlgIter, typename TWeight>
void extract_output(const TAlgIter& alg_it, 
        mxArray*& mxVs, mxArray*& mxDists, mxArray*& mxPreds)
{
    const std::vector<int>& vs = alg_it.close_sequence();
    int n = (int)vs.size();
    
    mxVs = iter_to_matlab_row(vs.begin(), n, 
            vnode_getter<TAlgIter>(alg_it));
    
    mxDists = iter_to_matlab_row(vs.begin(), n, 
            sdist_getter<TAlgIter, TWeight>(alg_it));
    
    mxPreds = iter_to_matlab_row(vs.begin(), n, 
            vpred_getter<TAlgIter>(alg_it));
}


template<typename T>
void do_solve(const MArray& mG, int s, char code, 
       mxArray*& mxVs, mxArray*& mxDists, mxArray*& mxPreds)
{
    RefWGraph<T> G = to_refwgraph<T>(mG);
    WAdjList<T> adjList = WAdjList<T>(G);

    if (code == 'd')
    {
        Dijkstra_SPathIterator<T> sp_it(adjList);
        sp_it.initialize(s);
        sp_it.solve_all();
                
        extract_output<Dijkstra_SPathIterator<T>, T>(sp_it, mxVs, mxDists, mxPreds);        
    }
    else if (code == 'a')
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

