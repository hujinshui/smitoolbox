/**********************************************************
 *
 *  topiclda_varinfer_cimp.cpp
 *
 *  C++ mex implementation of topiclda_varinfer
 *
 *  Created by Dahua Lin, on Feb 5, 2012
 *
 **********************************************************/

#include <bcslib/matlab/bcs_mex.h>
#include "specfuncs.h"

using namespace bcs;
using namespace bcs::matlab;


const int NUM_OBJ_ITEMS = 5;

struct DocObj
{
    double ell_theta;
    double ell_z;
    double ell_w;
    double ent_gamma;
    double ent_phi;
};


struct Doc
{
    Doc() 
    : nwords(0), words(0), counts(0) 
    {}
    
    int nwords;
    const int32_t *words;
    const double *counts;
};


inline double calc_sum_objv(const DocObj& s)
{
    return s.ell_theta + s.ell_z + s.ell_w + s.ent_gamma + s.ent_phi;
}


class Corpus
{
public:
    
    // I and J are zero-based indices
    Corpus(int n, int len, const int32_t *I, const int32_t *J, const double *C)
    : m_ndocs(n), m_docs(new Doc[n]), m_words(I), m_wcounts(C), m_max_count(0)
    {        
        // 1st-pass scan #words
        
        if (n > 1)
        {
            for (int i = 0; i < len; ++i)
            {
                ++ m_docs[J[i]].nwords;
            }
            
            int o = 0;
            
            for (int i = 0; i < n; ++i)
            {
                Doc& doc = m_docs[i];
                
                int nw = doc.nwords;
                if (nw > 0)
                {
                    doc.words = m_words + o;
                    doc.counts = m_wcounts + o;
                    o += nw;
                    
                    if (nw > m_max_count) 
                        m_max_count = nw;
                }
            }
        }
        else
        {
            m_docs->nwords = len;            
            m_docs->words = m_words;
            m_docs->counts = m_wcounts;
            
            m_max_count = m_docs->nwords;
        }                
    }    
    
    ~Corpus()
    {
        delete[] m_docs;
    }
    
    int ndocs() const { return m_ndocs; }
    
    int max_count() const { return m_max_count; }
    
    const Doc& doc(int i) const
    {
        return m_docs[i];
    }
            
private:
    int m_ndocs;
    Doc *m_docs;
    const int32_t *m_words;
    const double *m_wcounts;
    int m_max_count;
    
}; // end class Corpus



class Model
{
public:
    Model(int V, int K, const double *Beta, double alpha)
    : m_V(V), m_K(K), m_Beta(Beta), m_alpha(alpha)
    {
    }
    
    int V() const { return m_V; }
    
    int K() const { return m_K; }
    
    double beta(int v, int k) const { return m_Beta[v + m_V * k]; }
    
    double alpha() const { return m_alpha; }
    
private:
    int m_V;
    int m_K;
    const double *m_Beta;
    double m_alpha;
};


inline double calc_entropy(int n, const double *p)
{
    double e = 0;
    for (int i = 0; i < n; ++i) e += p[i] * std::log(p[i]);
    return -e;
}


void calc_objv(DocObj& obj, const Model& model, const Doc& doc, 
        const double *gamma, const double *Phi, double *temp)
{
    int K = model.K();
    double alpha = model.alpha();
    
    // pre-compute psi(gamma_k) - psi(gsum) --> temp
    
    double gsum = 0;
    for (int k = 0; k < K; ++k) gsum += gamma[k];
    
    double psi_gsum = smi::digamma(gsum);
    
    for (int k = 0; k < K; ++k) temp[k] = smi::digamma(gamma[k]) - psi_gsum;

    // ell_theta
    
    double a = 0;
    for (int k = 0; k < K; ++k) a += temp[k];
    
    obj.ell_theta = smi::gammaln(K * alpha) - K * smi::gammaln(alpha) +
            (alpha - 1) * a;
    
    
    // ell_z
    
    int nw = doc.nwords;
    
    obj.ell_z = 0;
    
    const double *phi = Phi;
    for (int i = 0; i < nw; ++i, phi += K)
    {
        a = 0;
        for (int k = 0; k < K; ++k)
        {
            a += phi[k] * temp[k];
        }
        obj.ell_z += doc.counts[i] * a;
    }    
    
    
    // ell_w
    
    obj.ell_w = 0;
    
    phi = Phi;
    for (int i = 0; i < nw; ++i, phi += K)
    {
        int v = doc.words[i];
        a = 0;
        for (int k = 0; k < K; ++k)
        {
            a += phi[k] * std::log(model.beta(v, k));
        }
        obj.ell_w += doc.counts[i] * a;
    }
    
    
    // ent_theta
               
    obj.ent_gamma = - smi::gammaln(gsum);
    for (int k = 0; k < K; ++k)
    {
        obj.ent_gamma += (smi::gammaln(gamma[k]) - (gamma[k] - 1) * temp[k]);
    }
    
    // ent_phi
    
    phi = Phi;
    obj.ent_phi = 0;
    for (int i = 0; i < nw; ++i, phi += K)
    {
        obj.ent_phi += doc.counts[i] * calc_entropy(K, phi);
    }
    
}



// Perform inference over a particular document
//
// return whether converged
//
// gamma need to be pre-initialized
//
// prev_gamma, exp_psi: at least K elements
//
bool infer_on_doc(const Model& model, const Doc& doc, 
        double *gamma, double *Phi, 
        int maxIter, double tol, 
        double *prev_gamma, double *exp_psi)
{
    int V = model.V();
    int K = model.K();
    double alpha = model.alpha();
    
    int nw = doc.nwords;
    bool converged = false;
        
    for (int t = 0; t < maxIter; ++t)
    {                        
        // store previous gamma values (for convergence test)
        for (int k = 0; k < K; ++k) 
        {
            exp_psi[k] = std::exp(smi::digamma(gamma[k]));
            prev_gamma[k] = gamma[k];
            gamma[k] = alpha;
        }
                
        double *phi = Phi;
        for (int i = 0; i < nw; ++i, phi += K)
        {                        
            int32_t v = doc.words[i];
            double c = doc.counts[i];
                                    
            // calculate phi 
            double sum_phi = 0;
            for (int k = 0; k < K; ++k)
            {                
                sum_phi += (phi[k] = model.beta(v, k) * exp_psi[k]);
            }
                        
            // normalize phi
            double nrm_coeff = 1.0 / sum_phi;            
            for (int k = 0; k < K; ++k) phi[k] *= nrm_coeff;
                                                            
            // accumulate phi to gamma
            for (int k = 0; k < K; ++k)
            {
                gamma[k] += c * phi[k];
            }                      
        }
                             
        // decide convergence (using L1-error)
        double err = 0;
        for (int k = 0; k < K; ++k)
        {
            err += std::abs(gamma[k] - prev_gamma[k]);
        }
        
        if (err < tol) 
        {
            converged = true;
            break;
        }
    }
           
    return converged;
}


inline void accum_phi(int K, double *sumPhi, const Doc& doc, double w, 
        const double *Phi)
{
    int nw = doc.nwords;
    for (int j = 0; j < nw; ++j, Phi += K)
    {
        int v = doc.words[j];
        double wc = w * doc.counts[j];
        double *sp = sumPhi + v * K;
        
        for (int k = 0; k < K; ++k)
        {
            sp[k] += wc * Phi[k];
        }
    }
}



void do_infer(const Model& model, const Corpus& corpus, const double *w,
        int maxIter, double tol, 
        double *Gamma, bool *converged, double *sumPhi, double *obj_vs)
{
    int K = model.K();
    double *gamma = Gamma;
    
    // allocate temporary memory
    
    double *Phi = new double[K * corpus.max_count()];
    double *prev_gamma = new double[K];
    double *exp_psi = new double[K];
    
    DocObj *objs = (DocObj*)obj_vs;
    
    for (int i = 0; i < corpus.ndocs(); ++i, gamma += K)
    {                
        const Doc& doc = corpus.doc(i);          
                       
        // per-document inference
                
        bool cvg = infer_on_doc(model, doc, gamma, 
                Phi, maxIter, tol, prev_gamma, exp_psi);
        converged[i] = cvg;
        
        // accumulate to sumPhi
        
        accum_phi(K, sumPhi, doc, w[i], Phi); 
        
        // calc objective
        
        double *temp = exp_psi;
        calc_objv(objs[i], model, doc, gamma, Phi, temp);
    }
    
    // release temporary memory
    
    delete[] Phi;
    delete[] prev_gamma;
    delete[] exp_psi;
}



/**
 * Inputs
 *   [0] Beta:          The word distributions of topics [V x K]
 *   [1] alpha:         The Dirichlet prior parameter [scalar]
 *   [2] w:             The document weights [vector of length n]
 *   [3] I:             The whole list of words [int32 zero-based]
 *   [4] J:             The whole list of documents [int32 zero-based]
 *   [5] C:             The whole list of word counts [double]
 *   [6] MaxIter:       The maximum number of iterations per document [double]
 *   [7] Tol:           Tolerance of change of gamma at convergence [double]
 *   [8] Gamma0:        The initial Gamma values [K x n]
 *
 * Outputs
 *   [0] Gamma:         The gamma vectors [K x n double]
 *   [1] converged:     Whether the variational inference converged [1 x n logical]
 *   [2] sumPhi:        The (weighted) sum of phi vectors [K x V double]
 *                      (for beta estimation)
 *   [3] ObjVs:         The itemized objective function
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mBeta(prhs[0]);
    const_marray mAlpha(prhs[1]);
    const_marray mW(prhs[2]);
    const_marray mI(prhs[3]);
    const_marray mJ(prhs[4]);
    const_marray mC(prhs[5]);
    const_marray mMaxIter(prhs[6]);
    const_marray mTol(prhs[7]);
    const mxArray *mxGamma0 = prhs[8];
    
    int V = mBeta.nrows();
    int K = mBeta.ncolumns();
    int n = (int)mxGetN(mxGamma0);
    int len = mI.nelems();
    
    const double *Beta = mBeta.data<double>();
    double alpha = mAlpha.get_scalar<double>();
    const double *w = mW.data<double>();
    
    const int32_t *I = mI.data<int32_t>();
    const int32_t *J = mJ.data<int32_t>();
    const double *C = mC.data<double>();
    
    int maxIter = (int)mMaxIter.get_scalar<double>();
    double tol = mTol.get_scalar<double>();
    
    // prepare output
            
    mxArray *mxGamma = mxDuplicateArray(mxGamma0);    
    marray mConverged = create_marray<bool>(1, n);
    marray mSumPhi = create_marray<double>(K, V);
    marray mObjVs = create_marray<double>(NUM_OBJ_ITEMS, n);
    
    double *Gamma = mxGetPr(mxGamma);
    bool *converged = mConverged.data<bool>();
    double *sumPhi = mSumPhi.data<double>();    
    double *objvs = mObjVs.data<double>();
    
    // main
    
    Model model(V, K, Beta, alpha);    
    Corpus corpus(n, len, I, J, C);
        
    do_infer(model, corpus, w, maxIter, tol, Gamma, converged, sumPhi, objvs);
    
    // output
    
    plhs[0] = mxGamma;
    plhs[1] = mConverged.mx_ptr();
    plhs[2] = mSumPhi.mx_ptr();   
    plhs[3] = mObjVs.mx_ptr();
}


BCSMEX_MAINDEF
        
        