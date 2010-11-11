/********************************************************************
 *
 *  intgroup.cpp
 *
 *  The C++ mex implementation of intgroup
 *
 *  Created by Dahua Lin, on June 6, 2010
 *
 ********************************************************************/

#include <mex.h>
#include <vector>

template<typename T>
void get_vs(const T *v, int ne, int& v0, int& v1)
{
    if (ne == 1)
    {
        v0 = 1;
        v1 = (int)(*v);
    }
    else
    {
        v0 = (int)(v[0]);
        v1 = (int)(v[1]);
    }
}


inline void get_range(const mxArray *mxRgn, int& v0, int& v1)
{
    int ne = mxGetNumberOfElements(mxRgn);
    
    if (!( (ne == 1 || ne == 2) && !mxIsSparse(mxRgn)))
    {
        mexErrMsgIdAndTxt("intcount:invalidarg", 
                "The range [v0 v1] should be a (non-sparse) pair.");
    }
    
    if (mxIsDouble(mxRgn))
    {
        const double *v = (const double*)mxGetData(mxRgn);
        get_vs(v, ne, v0, v1);
    }
    else if (mxIsSingle(mxRgn))
    {
        const float *v = (const float*)mxGetData(mxRgn);
        get_vs(v, ne, v0, v1);
    }
    else if (mxIsInt32(mxRgn))
    {
        const int *v = (const int*)mxGetData(mxRgn);
        get_vs(v, ne, v0, v1);
    }
    else
    {
        mexErrMsgIdAndTxt("intcount:invalidarg", 
                "The range [v0 v1] should be of class double, single, or int32");
    }
}


template<typename T>
void group(int v0, int v1, const T *v, int n, std::vector<double>* gs)
{
    for (int i = 0; i < n; ++i)
    {
        int cv = (int)(v[i]);
        if (cv >= v0 && cv <= v1)
        {
            gs[cv - v0].push_back(i + 1);
        }
    }
}

inline void do_group(int v0, int v1, const mxArray *mxVals, std::vector<double>* gs)
{        
    
    int n = mxGetNumberOfElements(mxVals);
    
    if (mxIsDouble(mxVals))
    {
        const double *v = (const double*)mxGetData(mxVals);
        group(v0, v1, v, n, gs);
    }
    else if (mxIsSingle(mxVals))
    {
        const float *v = (const float*)mxGetData(mxVals);
        group(v0, v1, v, n, gs);
    }
    else if (mxIsInt32(mxVals))
    {
        const int *v = (const int*)mxGetData(mxVals);
        group(v0, v1, v, n, gs);
    }
    else
    {
        mexErrMsgIdAndTxt("intgroup:invalidarg", 
                "The class of values should be either double, single, int32");
    }    
}


mxArray *to_mx_vector(const std::vector<double>& g)
{
    int n = (int)g.size();
    mxArray *mxV = mxCreateDoubleMatrix(1, n, mxREAL);
    
    double *v = mxGetPr(mxV);
    for (int i = 0; i < n; ++i)
    {
        v[i] = g[i];
    }
    
    return mxV;
}


mxArray *groups_to_cells(int ng, const std::vector<double>* gs)
{    
    mxArray *mxC = mxCreateCellMatrix(1, ng);
    
    for (int i = 0; i < ng; ++i)
    {
        mxSetCell(mxC, i, to_mx_vector(gs[i]));
    }
    
    return mxC;    
}


/**
 * Main entry:
 *
 *   Input:
 *     [0]: rgn:    the integer range in form of [v0, v1] (pair) 
 *     [1]: V:      the values array (double|single|int32)
 *   Output:
 *     [0]: G:      the cell array of indices for each value
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
        mexErrMsgIdAndTxt("intgroup:invalidarg", 
                "The number of inputs to intgroup should be 2.");
    
    const mxArray *mxRgn = prhs[0];
    const mxArray *mxVals = prhs[1];
    
    int v0, v1;
    get_range(mxRgn, v0, v1);    
    
    int ng = v1 - v0 + 1;
    std::vector<double> *vecs = new std::vector<double>[ng];
    
    do_group(v0, v1, mxVals, vecs);            
    
    plhs[0] = groups_to_cells(ng, vecs);
    
    delete[] vecs;
}




