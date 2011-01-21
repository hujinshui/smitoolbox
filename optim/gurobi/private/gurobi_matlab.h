/********************************************************************
 *
 *  gurobi_matlab.h
 *
 *  A light-weight wrapper of Gurobi model for mex
 *
 *  Created by Dahua Lin, on Jan 20, 2011
 *
 ********************************************************************/

#include <mex.h>

extern "C"
{
#include <gurobi_c.h>
}

#include <string.h>
#include <vector>

namespace gurobi
{
    
// an exception class that captures Gurobi error    
class Error
{
public:
    explicit Error(int code) : _code(code)
    {
    }
    
    int code() const
    {
        return _code;
    }
    
    const char *message() const
    {
        switch (_code)
        {
            case 0:
                return "No error.";
            case GRB_ERROR_OUT_OF_MEMORY:
                return "Available memory was exhausted.";
            case GRB_ERROR_NULL_ARGUMENT:
                return "NULL input value was provided for a required argument.";
            case GRB_ERROR_INVALID_ARGUMENT:
                return "An invalid value was provided for a routine argument.";
            case GRB_ERROR_UNKNOWN_ATTRIBUTE:
                return "Tried to query or set an unknown attribute.";
            case GRB_ERROR_DATA_NOT_AVAILABLE:
                return "Attempted to query or set an attribute that could not be accessed.";                
            case GRB_ERROR_INDEX_OUT_OF_RANGE:
                return "One or more of the provided indices was out of range.";
            case GRB_ERROR_UNKNOWN_PARAMETER:
                return "Tried to query or set an unknown parameters.";                
            case GRB_ERROR_VALUE_OUT_OF_RANGE:
                return "Tried to set a value that is out of range.";
            case GRB_ERROR_NO_LICENSE:
                return "Failed to obtain a valid Gurobi license.";
            case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
                return "Attempted to solve a model larger than the size allowed in license.";
            case GRB_ERROR_CALLBACK:
                return "Encountered problem in callback.";
            case GRB_ERROR_FILE_READ:
                return "Failed to read an requested file.";
            case GRB_ERROR_FILE_WRITE:
                return "Failed to write an requested file.";
            case GRB_ERROR_NUMERIC:
                return "Encountered numeric error during an operation";
            case GRB_ERROR_IIS_NOT_INFEASIBLE:
                return "Attempted to perform infeasibility analysis in a feasible model";
            case GRB_ERROR_NOT_FOR_MIP:
                return "The requested operation is not valid for an MIP model.";
            case GRB_ERROR_OPTIMIZATION_IN_PROGRESS:
                return "Tried to query or modify a model when the optimization was in progress.";
            case GRB_ERROR_DUPLICATES:
                return "Constraints, variables, or SOS contained duplicated indices.";
            case GRB_ERROR_NODEFILE:
                return "Error in reading or writing a node file during MIP optimization.";
            case GRB_ERROR_Q_NOT_PSD:
                return "The quadratic matrix in QP model is not positive semi-definite";
            default:
                return "Unknown error.";                
        }
    }
    
private:
    int _code;
};
    
    
    
// a light-weight wrapper of Gurobi environment    
class Env
{
public:
    explicit Env() : _env(0), _standalone(true)
    {
        int err = GRBloadenv(&_env, NULL);
        if (err != 0)
            throw Error(err);            
    }
    
    Env(GRBenv *e) : _env(e), _standalone(false)
    {
    }
    
    ~Env()
    {
        if (_standalone)
        {
            GRBfreeenv(_env);
        }
        
        if (!_strs.empty())
        {
            int n = (int)_strs.size();
            for (int i = 0; i < n; ++i)
            {
                delete[] _strs[i];
            }
        }
    }
    
    GRBenv *ptr()
    {
        return _env;
    }
    
    void set_int_param(const char *name, int v)
    {
        int err = GRBsetintparam(_env, name, v);
        if (err != 0)
            throw Error(err);
    }
    
    void set_double_param(const char *name, double v)
    {        
        int err = GRBsetdblparam(_env, name, v);
        if (err != 0)
            throw Error(err);
    }
    
    void set_string_param(const char *name, const char *v)
    {
        int len = ::strlen(v);
        char *sz = new char[len + 1];
        ::memcpy(sz, v, (len+1) * sizeof(char));
        _strs.push_back(sz);
        
        int err = GRBsetstrparam(_env, name, sz);
        if (err != 0)
            throw Error(err);                
    }
    
private:
    Env(const Env& );
    Env& operator = (const Env& );
    
private:
    bool _standalone;
    GRBenv *_env;
    std::vector<char*> _strs;
};
    
    
    
// a light-weight wrapper of Gurobi model    
class Model
{
public:
    explicit Model(Env& env, const char *name = NULL) 
    : _env(env), _model(0), _nvars(0), _ncons(0)
    {
        int err = GRBnewmodel(env.ptr(), &_model, name, 0, 
                NULL, NULL, NULL, NULL, NULL);
        if (err != 0)
            throw Error(err);
    }
    
    ~Model()
    {
        GRBfreemodel(_model);
        
        if (!_rels.empty())
        {
            int nr = (int)_rels.size();
            for (int i = 0; i < nr; ++i)
            {
                char *pr = _rels[i];
                if (pr) delete[] pr;                
            }
            _rels.clear();
        }
    }
    
    
    int numvars() const
    {
        return _nvars;
    }
    
    int numconstrs() const
    {
        return _ncons;
    }
    
    GRBmodel *ptr() 
    {
        return _model;
    }
    
    
    
    /**
     * add variables to the model
     *
     * @param n  the number of variables to be added
     * @param c  the coefficients in objective
     * @param lb the lower bounds of the variables (can be NULL)
     * @param ub the upper bounds of the variables (can be NULL)
     *     
     */
    void addvars(int n, 
            const double *c, const double *lb, const double *ub)
    {
        double *c_ = const_cast<double*>(c);
        double *lb_ = const_cast<double*>(lb);
        double *ub_ = const_cast<double*>(ub);
                        
        int err = GRBaddvars(_model, n, 0, NULL, NULL, NULL, 
                c_, lb_, ub_, NULL, NULL);
        
        if (err == 0)
        {
            _nvars += n;
        }
        else
        {
            throw Error(err);
        }                    
    }
    
    
    /**
     * add constraints to the model
     *
     * @param m   the number of constraints to be added
     * @param nnz the number of non-zero coefficients in these constraints
     * @param rel the relation {GRB_LESS_EQUAL, GRB_GREATER_EQUAL, GRB_EQUAL}
     * @param os  the offsets for each constraint
     * @param ind the column indices of non-zero coefficients for each constraint
     * @param c   the non-zero coefficients for each constraint
     * @param rhvs the right hand side values
     *
     */
    void addconstrs(int m, int nnz, char rel, 
            const int *os, const int *ind, const double *c, const double *rhvs)
    {
        int pm = _ncons;
     
        int *cbeg = const_cast<int*>(os);
        int *cind = const_cast<int*>(ind);
        double *cval = const_cast<double*>(c);
                
        char *sense = new char[m];
        for (int i = 0; i < m; ++i) sense[i] = rel;
        
        double *rhs = const_cast<double*>(rhvs);
        
        int err = GRBaddconstrs(_model, m, nnz, cbeg, cind, cval, 
                sense, rhs, NULL);
        
        if (err == 0)
        {
            _ncons += m;
            _rels.push_back(sense);
        }
        else
        {
            delete[] sense;
            throw Error(err);
        }
    }
    
    
    /**
     * add quadratic terms for QP
     *
     * @param n     the number of terms
     * @param I     the array of row indices (zero-based)
     * @param J     the array of column indices (zero-based)
     * @param V     the coefficient values
     */
    void addqpterms(int n, const int *I, const int* J, const double *V)
    {
        int err = GRBaddqpterms(_model, n,
                const_cast<int*>(I),
                const_cast<int*>(J),
                const_cast<double*>(V));
        
        if (err != 0)
            throw Error(err);
    }
    
    
    /**
     * Update the model to integrate previous modifications
     */
    void update()
    {
        int err = GRBupdatemodel(_model);
        
        if (err != 0)
            throw Error(err);
    }
    
    
    /**
     * Do optimization
     */
    void solve()
    {
        int err = GRBoptimize(_model);
        
        if (err != 0)
            throw Error(err);
    }        
    
    
    void set_callback(int (*cb)(GRBmodel*,void *,int,void*), void *usrdata=NULL)
    {
        GRBsetcallbackfunc(_model, cb, usrdata);
    }
    
    
    int get_solution(double *x) const
    {
        GRBmodel *model = const_cast<GRBmodel*>(_model);        
        return GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, _nvars, x);
    }
    
    
    /******************************************
     *
     *  Attribute getting
     *
     ******************************************/
    
    
    int query_objval(double &v) const
    {
        return query_double_attr(GRB_DBL_ATTR_OBJVAL, v);
    }
    
    int query_objbound(double &v) const
    {
        return query_double_attr(GRB_DBL_ATTR_OBJBOUND, v);
    }
    
    int query_status(int &v) const
    {
        return query_int_attr(GRB_INT_ATTR_STATUS, v);
    }
    
    int query_runtime(double &v) const
    {
        return query_double_attr(GRB_DBL_ATTR_RUNTIME, v);
    }
    
    int query_itercount(double& v) const
    {
        return query_double_attr(GRB_DBL_ATTR_ITERCOUNT, v);
    }
    
    int query_baritercount(int& v) const
    {
        return query_int_attr(GRB_INT_ATTR_BARITERCOUNT, v);
    }
    
    
    
    int query_int_attr(const char *name, int& v) const
    {
        GRBmodel *model = const_cast<GRBmodel*>(_model);        
        return GRBgetintattr(model, name, &v);
    }
    
    int query_double_attr(const char *name, double& v) const
    {
        GRBmodel *model = const_cast<GRBmodel*>(_model);
        return GRBgetdblattr(model, name, &v);
    }            
    
private:
    Env& _env;
    GRBmodel *_model;
    
    int _nvars;     // the number of variables
    int _ncons;     // the number of constraints
    std::vector<char*> _rels;    // the constraint relations
    
}; // end class Model
    
    
    
} // end namespace gurobi


