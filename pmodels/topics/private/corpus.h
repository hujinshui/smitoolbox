// The header file for the class Corpus

#ifndef SMI_CORPUS_H
#define SMI_CORPUS_H


namespace smi
{


struct Doc
{
    Doc() 
    : nwords(0), words(0), counts(0) 
    {}
    
    int nwords;
    const int32_t *words;
    const double *counts;
};


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


} // end namespace smi


#endif
