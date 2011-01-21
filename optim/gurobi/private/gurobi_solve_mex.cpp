/********************************************************************
 *
 *  gurobi_solve_mex.cpp
 *
 *  The mex wrapper of Gurobi solver
 *
 *  Created by Dahua Lin, on Jan 21, 2011
 *
 ********************************************************************/

#include "gurobi_matlab.h"

using namespace gurobi;

/**
 * The class to represent a constraint matrix (sparse)
 */
class ConsMat
{
public:
    
    ConsMat()
    : _m(0), _N(0), _cbeg(0), _cind(0), _cval(0)
    {
    }
    
    ~ConsMat()
    {
        if (_cbeg != 0)
        {
            delete[] _cbeg;
            delete[] _cind;
            delete[] _cval;
        }        
    }
    
    bool empty() const
    {
        return _m == 0;
    }
    
                    
    int nc() const
    {
        return _m;
    }
    
    int nnz() const
    {
        return _N;
    }
    
        
    const int *cbeg() const
    {
        return _cbeg;
    }
    
    const int *cind() const
    {
        return _cind;
    }
    
    const double *cval() const
    {
        return _cval;
    }
    
    void create_from(const mxArray *mxS)  // can be invoked only once 
    {
        const mxArray *mxM = mxGetField(mxS, 0, "m");
        const mxArray *mxI = mxGetField(mxS, 0, "i");
        const mxArray *mxJ = mxGetField(mxS, 0, "j");
        const mxArray *mxV = mxGetField(mxS, 0, "v");
        
        _m = (int)mxGetScalar(mxM);
        _N = mxGetNumberOfElements(mxV);
        
        const int *I = (const int*)mxGetData(mxI);
        const int *J = (const int*)mxGetData(mxJ);
        const double *V = (const double*)mxGetData(mxV);
        
        _cbeg = new int[_m];
        _cind = new int[_N];
        _cval = new double[_N];
        
        // 1st scan: fill cbeg
        
        for (int k = 0; k < _m; ++k) _cbeg[k] = 0;
        for (int i = 0; i < _N; ++i) ++ _cbeg[I[i]];        
        
        int b = 0;
        for (int k = 0; k < _m; ++k) 
        {
            int b_p = b;
            b += _cbeg[k];
            _cbeg[k] = b_p;
        }
        
        // 2nd scan: fill cind and cval
        int *os = new int[_m];
        for (int k = 0; k < _m; ++k) os[k] = 0;
        
        for (int i = 0; i < _N; ++i)
        {
            int k = I[i];
            int j = J[i];
            double v = V[i];
            
            int o = _cbeg[k] + os[k];
            
            _cind[o] = j;
            _cval[o] = v;
            
            ++ os[k];
        }
                
        delete[] os;
    }
    
    
private:
    int _m;
    int _N;
    
    int *_cbeg;
    int *_cind;
    double *_cval;
    
}; // end class ConsMat



/**
 * The class to represent a quadratic matrix
 */
class QMat
{
public:
    QMat(const mxArray *mxQ)
    : _nnz(0), _I(0), _J(0), _V(0)
    {        
        const mxArray *mxI = mxGetField(mxQ, 0, "I");
        const mxArray *mxJ = mxGetField(mxQ, 0, "J");
        const mxArray *mxV = mxGetField(mxQ, 0, "V");
        
        _nnz = mxGetNumberOfElements(mxI);
        
        _I = (const int*)mxGetData(mxI);
        _J = (const int*)mxGetData(mxJ);
        _V = mxGetPr(mxV);
    }
    
    int nnz() const
    {
        return _nnz;
    }
    
    const int *qrows() const
    {
        return _I;
    }
    
    const int *qcols() const
    {
        return _J;
    }
    
    const double *qvals() const
    {
        return _V;
    }
            
private:
    int _nnz;
    const int *_I;
    const int *_J;
    const double *_V;
    
}; // end class QMat



inline mxArray *empty_mat()
{
    return mxCreateDoubleMatrix(0, 0, mxREAL);
}


inline int status_to_flag(int s)
{
    if (s == GRB_OPTIMAL)
    {
        return 0;
    }
    else if (s == GRB_INF_OR_UNBD || s == GRB_INFEASIBLE || s == GRB_UNBOUNDED)
    {
        return 1;
    }
    else
    {
        return 2;
    }
}


struct AddOptions
{
public:
    AddOptions()
    {
        to_display = false;
    }
        
    bool to_display;      
};


void set_params(Env& env, const mxArray *mxParams, AddOptions& oa)
{            
    int nf = (int)mxGetNumberOfFields(mxParams);
    for (int i = 0; i < nf; ++i)
    {
        const char *name = mxGetFieldNameByNumber(mxParams, i);                
        const mxArray *mxV = mxGetField(mxParams, 0, name);
                                
        if (mxIsInt32(mxV))
        {            
            int v = *((const int*)mxGetData(mxV));
            
            if (::strcmp(name, "Display") == 0)
            {
                oa.to_display = (v != 0);
            }
            else
            {
                env.set_int_param(name, v); 
            }                                   
        }
        else if (mxIsDouble(mxV))
        {
            double v = mxGetScalar(mxV);
            env.set_double_param(name, v);
        }
        else if (mxIsChar(mxV))
        {
            int len = mxGetNumberOfElements(mxV) + 1;
            char *sz = new char[len];
            mxGetString(mxV, sz, len);
            sz[len] = '\0';
            
            env.set_string_param(name, sz);
            
            delete[] sz;
        }         
    }
     
}



mxArray *make_info_struct(int status, const Model& model)
{
    static const char* fieldnames[] = {
        "status", "runtime", "iters", "bariters", "objval", "objbound"
    };    
    static int nfields = sizeof(fieldnames) / sizeof(const char*);
    
    double runtime;
    double iters;
    int bariters;
    double objval;
    double objbound;   
        
    mxArray *mxS = mxCreateStructMatrix(1, 1, nfields, fieldnames);
        
    mxSetField(mxS, 0, "status", mxCreateDoubleScalar(status));
    
    if (model.query_runtime(runtime) == 0)
    {
        mxSetField(mxS, 0, "runtime", mxCreateDoubleScalar(runtime));
    }
    if (model.query_itercount(iters) == 0)
    {
        mxSetField(mxS, 0, "iters", mxCreateDoubleScalar(iters));
    }
    if (model.query_baritercount(bariters) == 0)
    {
        mxSetField(mxS, 0, "bariters", mxCreateDoubleScalar(bariters));
    }
    if (model.query_objval(objval) == 0)
    {
        mxSetField(mxS, 0, "objval",   mxCreateDoubleScalar(objval));
    }
    if (model.query_objbound(objbound) == 0)
    {
        mxSetField(mxS, 0, "objbound", mxCreateDoubleScalar(objbound));
    }                    
        
    return mxS;
}



// message handling callback function
int __stdcall my_callback(GRBmodel *model, void *cbdata, int where, void *usrdata)
{
    char *msg;
    
    if (where == GRB_CB_MESSAGE)
    {
        int err = GRBcbget(cbdata, where, GRB_CB_MSG_STRING, &msg);
        if (err != 0)
            throw Error(err);
        mexPrintf("%s", msg);
    }    
        
    return 0;
}




// main entry
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take input
    
    const mxArray *mxS = prhs[0];
    const mxArray *mxParams = prhs[1];
    
    const mxArray *mxD  = mxGetField(mxS, 0, "d");
    const mxArray *mxF  = mxGetField(mxS, 0, "f");
    
    const mxArray *mxA  = mxGetField(mxS, 0, "A");
    const mxArray *mxBl = mxGetField(mxS, 0, "bl");
    const mxArray *mxBu = mxGetField(mxS, 0, "bu");
    const mxArray *mxAe = mxGetField(mxS, 0, "Ae");
    const mxArray *mxBe = mxGetField(mxS, 0, "be");
    
    const mxArray *mxLb = mxGetField(mxS, 0, "lb");
    const mxArray *mxUb = mxGetField(mxS, 0, "ub");
    
    const mxArray *mxQ = mxGetField(mxS, 0, "Q");
    
    // process input
    
    int n = (int)mxGetScalar(mxD);
    const double *f = mxGetPr(mxF);
    
    ConsMat A;    
    const double *bl = 0;
    const double *bu = 0;    
    
    if (!mxIsEmpty(mxA))
    {
        A.create_from(mxA);
        bl = mxGetPr(mxBl);
        bu = mxGetPr(mxBu);
    }
    
    ConsMat Ae;
    const double *be = 0;
    
    if (!mxIsEmpty(mxAe))
    {
        Ae.create_from(mxAe);
        be = mxGetPr(mxBe);
    }
    
    const double *lb = 0;
    const double *ub = 0;
    
    if (!mxIsEmpty(mxLb))
    {
        lb = mxGetPr(mxLb);
    }
    
    if (!mxIsEmpty(mxUb))
    {
        ub = mxGetPr(mxUb);
    }
    
    
    // main

    try
    {
        // prepare model
        
        Env env;               
        Model model(env);
        
        model.addvars(n, f, lb, ub);        
        model.update();
        
        if (!mxIsEmpty(mxQ))
        {
            QMat Q(mxQ);
            model.addqpterms(Q.nnz(), Q.qrows(), Q.qcols(), Q.qvals());
        }        
        
        if (bl != 0)
        {
            model.addconstrs(A.nc(), A.nnz(), GRB_GREATER_EQUAL, 
                    A.cbeg(), A.cind(), A.cval(), bl);
        }
        
        if (bu != 0)
        {
            model.addconstrs(A.nc(), A.nnz(), GRB_LESS_EQUAL,
                    A.cbeg(), A.cind(), A.cval(), bu);
        }
        
        if (be != 0)
        {
            model.addconstrs(Ae.nc(), Ae.nnz(), GRB_EQUAL,
                    Ae.cbeg(), Ae.cind(), Ae.cval(), be);
        }
        
        model.update();
        
        // set parameters
        
        AddOptions oa;
        
        if (!mxIsEmpty(mxParams))
        {
            Env env_a(GRBgetenv(model.ptr()));                
            set_params(env_a, mxParams, oa);
        }                                
        
        // solve
        
        if (oa.to_display)
        {
            model.set_callback(my_callback);
        }
        
        model.solve();
        
        // get solution
        
        int status = 0;
        model.query_status(status);
        
        if (status == GRB_OPTIMAL)
        {
            mxArray *mxSol = mxCreateDoubleMatrix(n, 1, mxREAL);
            model.get_solution(mxGetPr(mxSol));
            plhs[0] = mxSol;
        }
        else
        {
            plhs[0] = empty_mat();
        }
        
        // get objective value
        
        if (nlhs >= 2)
        {
            if (status == GRB_OPTIMAL)
            {
                double objv = 0;
                model.query_objval(objv);
                plhs[1] = mxCreateDoubleScalar(objv);
            }
            else
            {
                plhs[1] = empty_mat();
            }
        }
        
        // get flag
        
        if (nlhs >= 3)
        {
            int flag = status_to_flag(status);            
            plhs[2] = mxCreateDoubleScalar(flag);            
        }
        
        // get info
        
        if (nlhs >= 4)
        {
            plhs[3] = make_info_struct(status, model);
        }                        
        
    }
    catch(Error err)
    {
        mexErrMsgIdAndTxt("gurobi_solve:exception", err.message());
    } 
    
}












