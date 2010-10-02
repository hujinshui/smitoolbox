/********************************************************************
 *
 *  gr_neighbors_cimp.cpp
 *
 *  The C++ mex implementation of gr_neighbor.m
 *
 *  Created by Dahua Lin, on Oct 1, 2010
 *
 ********************************************************************/

#include "../../../base/clib/matlab_types.h"
#include "../../clib/mgraph.h"

using namespace smi;

template class MWGraph<double>;
template class WGNeighborHood<double>;

inline void raise_operr()
{
    mexErrMsgIdAndTxt("gr_neighbors:invalidarg", "The op char is invalid.");
}


mxArray* nbnodes_to_matlab_cells(const GNeighborHood& nbh)
{
    int n = nbh.nnodes();
    mxArray *mxC = mxCreateCellMatrix(n, 1);
    for (int i = 0; i < n; ++i)
    {
        mxArray *mxI = src_to_matlab_matrix<int>(1, nbh.neighbor_num(i), nbh.neighbor_nodes(i));
        mxSetCell(mxC, i, mxI);
    }    
    return mxC;
}

template<typename T>
mxArray* nbweights_to_matlab_cells(const WGNeighborHood<T>& nbh)
{
    int n = nbh.nnodes();
    mxArray *mxC = mxCreateCellMatrix(n, 1);
    for (int i = 0; i < n; ++i)
    {
        mxArray *mxI = src_to_matlab_matrix<T>(1, nbh.neighbor_num(i), nbh.neighbor_weights(i));
        mxSetCell(mxC, i, mxI);
    }
    return mxC;
}


mxArray* do_extract_nbs(const mxArray *mxG, char op)
{
    MGraph G(mxG);
    
    if (op == 'o' || op == 'O')
    {
        GNeighborHood nbh(G, gnb_out());
        return nbnodes_to_matlab_cells(nbh);        
    }
    else if (op == 'i' || op == 'I')
    {
        GNeighborHood nbh(G, gnb_in());
        return nbnodes_to_matlab_cells(nbh);
    }
    else
    {
        raise_operr();
        return 0;
    }
}


template<typename T>
void do_extract_wnbs(const mxArray *mxG, char op, mxArray*& mxNbs, mxArray*& mxWs)
{
    MWGraph<T> G(mxG);
    
    if (op == 'o' || op == 'O')
    {
        WGNeighborHood<T> nbh(G, gnb_out());
        mxNbs = nbnodes_to_matlab_cells(nbh);     
        mxWs = nbweights_to_matlab_cells(nbh);
    }
    else if (op == 'i' || op == 'I')
    {
        WGNeighborHood<T> nbh(G, gnb_in());
        mxNbs = nbnodes_to_matlab_cells(nbh);
        mxWs = nbweights_to_matlab_cells(nbh);
    }
    else
    {
        raise_operr();
    }
}
       

// main entry:
//  Input:
//    [0] G:    the mgraph struct
//    [1] op:   the option code (char)
//
//  Output:
//    [0] nbs:  the cell array of neighbor indices
//    [1] nws:  the cell array of neighbor weights (optional)
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *mxG = prhs[0];
    const mxArray *mxOp = prhs[1];
    
    char op = (char)(*((const mxChar*)mxGetData(mxOp)));
    
    if (nlhs <= 1)
    {
        plhs[0] = do_extract_nbs(mxG, op);
    }
    else
    {
        mxClassID wcid = MGraph::get_weight_class(mxG);
        
        switch (wcid)
        {
            case mxDOUBLE_CLASS:
                do_extract_wnbs<double>(mxG, op, plhs[0], plhs[1]);
                break;
            case mxINT32_CLASS:
                do_extract_wnbs<int>(mxG, op, plhs[0], plhs[1]);
                break;
            case mxSINGLE_CLASS:
                do_extract_wnbs<float>(mxG, op, plhs[0], plhs[1]);
                break;
        }        
    }            
    
}


