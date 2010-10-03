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

#include <mex.h>
#include "graphs.h"

namespace smi
{
    
inline mxClassID get_graph_weight_class(const mxArray *mxG)
{
    const mxArray *mxW = mxGetField(mxG, 0, "W");
    return mxGetClassID(mxW);
}
    
    
inline RefGraph to_refgraph(const mxArray *mxG)
{
    const mxArray *mxN = mxGetField(mxG, 0, "n");
    const mxArray *mxI = mxGetField(mxG, 0, "I");
    const mxArray *mxJ = mxGetField(mxG, 0, "J");        
    
    int n = (int)mxGetScalar(mxN);
    int m = mxGetNumberOfElements(mxI);
    const int *I = (const int*)mxGetData(mxI);
    const int *J = (const int*)mxGetData(mxJ);
    
    return RefGraph(n, m, I, J);    
}


template<typename TWeight>
inline RefWGraph<TWeight> to_refwgraph(const mxArray *mxG)
{
    const mxArray *mxN = mxGetField(mxG, 0, "n");
    const mxArray *mxI = mxGetField(mxG, 0, "I");
    const mxArray *mxJ = mxGetField(mxG, 0, "J"); 
    const mxArray *mxW = mxGetField(mxG, 0, "W");
    
    int n = (int)mxGetScalar(mxN);
    int m = mxGetNumberOfElements(mxI);
    const int *I = (const int*)mxGetData(mxI);
    const int *J = (const int*)mxGetData(mxJ);
    const TWeight *W = (const TWeight*)mxGetData(mxW);
    
    return RefWGraph<TWeight>(n, m, I, J, W);
}
        

inline mxArray *gindices_mrow(int n, int *v)
{
    mxArray *mxI = mxCreateNumericMatrix(1, n, mxINT32_CLASS, mxREAL);
    int *I = (int*)mxGetData(mxI);
    for (int i = 0; i < n; ++i)
    {
        I[i] = v[i] + 1;
    }
    return mxI;
}

}

#endif





