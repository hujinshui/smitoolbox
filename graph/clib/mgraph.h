/********************************************************************
 *
 *  mgraph.h
 *
 *  The interfaces between matlab and C++ classes
 *
 *  Created by Dahua Lin, on Oct 3, 2010
 *
 *******************************************************************/

#ifndef SMI_MGRAPH_H
#define SMI_MGRAPH_H

#include "../../base/clib/marray.h"
#include "graphs.h"

namespace smi
{
    
inline mxClassID get_graph_weight_class(const MArray& mG)
{
    return mG.get_field("W").class_id();   
}
    
    
inline RefGraph to_refgraph(const MArray& mG)
{
    int n = (int)mG.get_field("n").get_double_scalar();
    MArray mI = mG.get_field("I");
    MArray mJ = mG.get_field("J");        
    int m = mI.nelems();
    
    const int *I = mI.get_data<int>();
    const int *J = mJ.get_data<int>();
    
    return RefGraph(n, m, I, J);    
}


template<typename TWeight>
inline RefWGraph<TWeight> to_refwgraph(const MArray &mG)
{
    int n = (int)mG.get_field("n").get_double_scalar();
    MArray mI = mG.get_field("I");
    MArray mJ = mG.get_field("J"); 
    MArray mW = mG.get_field("W");
    int m = mI.nelems();
    
    const int *I = mI.get_data<int>();
    const int *J = mJ.get_data<int>();
    const TWeight *W = mW.get_data<TWeight>();
    
    return RefWGraph<TWeight>(n, m, I, J, W);
}
        

inline mxArray *gindices_mrow(int n, const int *v)
{
    mxArray *mxI = create_matlab_matrix<int>(1, n);
    int *I = (int*)mxGetData(mxI);
    for (int i = 0; i < n; ++i)
    {
        I[i] = v[i] + 1;
    }
    return mxI;
}

}

#endif





