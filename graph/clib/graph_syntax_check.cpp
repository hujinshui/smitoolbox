// A special file for syntax checking of the template classes

#include <mex.h>

#include "graph_base.h"
#include "graph_adjlist.h"

using namespace smi;


template <class G>
struct RefVertexListGraphConcept 
{
    typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
    
    void constraints() 
    {
        G& g = *pg;
        
        boost::function_requires< boost::MultiPassInputIteratorConcept<vertex_iterator> >();
        
        p = vertices(g);
        V = num_vertices(g);
        v = *p.first;
        const_constraints(g);
    }
    
    void const_constraints(const G& g) 
    {
        p = vertices(g);
        V = num_vertices(g);
        v = *p.first;
    }
    
    std::pair<vertex_iterator, vertex_iterator> p;
    typename boost::graph_traits<G>::vertex_descriptor v;
    typename boost::graph_traits<G>::vertices_size_type V;
    G *pg;
};


template <class G>
struct RefEdgeListGraphConcept 
{
    typedef typename boost::graph_traits<G>::edge_iterator edge_iterator;
    
    void constraints() 
    {
        boost::function_requires< boost::GraphConcept<G> >();
        boost::function_requires< boost::MultiPassInputIteratorConcept<edge_iterator> >();
        
        G& g = *pg;
        
        p = edges(g);
        E = num_edges(g);
        e = *p.first;
        u = source(e, g);
        v = target(e, g);
        const_constraints(g);
    }
    
    void const_constraints(const G& g) 
    {
        p = edges(g);
        E = num_edges(g);
        e = *p.first;
        u = source(e, g);
        v = target(e, g);
    }
    
    std::pair<edge_iterator, edge_iterator> p;
    typename boost::graph_traits<G>::vertex_descriptor u, v;
    typename boost::graph_traits<G>::edge_descriptor e;
    typename boost::graph_traits<G>::edges_size_type E;
    G *pg;
};


template <class G>
struct RefIncidenceGraphConcept 
{
    typedef typename boost::graph_traits<G>::out_edge_iterator out_edge_iterator;
    
    void constraints() 
    {
        boost::function_requires< boost::GraphConcept<G> >();
        boost::function_requires< boost::MultiPassInputIteratorConcept<out_edge_iterator> >();
        
        G& g = *pg;
        
        p = out_edges(u, g);
        e = *p.first;
        u = source(e, g);
        v = target(e, g);
    }
    
    void const_constraints(const G& g) 
    {
        p = out_edges(u, g);
        e = *p.first;
        u = source(e, g);
        v = target(e, g);
    }
    
    std::pair<out_edge_iterator, out_edge_iterator> p;
    typename boost::graph_traits<G>::vertex_descriptor u, v;
    typename boost::graph_traits<G>::edge_descriptor e;
    G *pg;
};


template <class G>
struct RefAdjacencyGraphConcept 
{
    typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
    
    void constraints() 
    {
        boost::function_requires< boost::GraphConcept<G> >();
        boost::function_requires< boost::MultiPassInputIteratorConcept<adjacency_iterator> >();
        
        G& g = *pg;
        p = adjacent_vertices(v, g);
        v = *p.first;
        const_constraints(g);
    }
    
    void const_constraints(const G& g) 
    {
        p = adjacent_vertices(v, g);
    }
    
    std::pair<adjacency_iterator, adjacency_iterator> p;
    typename boost::graph_traits<G>::vertex_descriptor v;
    G *pg;
};






template<typename TWeight>
struct CRefEdgeListConceptCheck
{
    BOOST_CONCEPT_ASSERT((RefVertexListGraphConcept<CRefEdgeList<TWeight> >)); 
    BOOST_CONCEPT_ASSERT((RefEdgeListGraphConcept<CRefEdgeList<TWeight> >));   
};

template<typename TWeight>
struct CRefAdjListConceptCheck
{
    BOOST_CONCEPT_ASSERT((RefVertexListGraphConcept<CRefAdjList<TWeight> >)); 
    BOOST_CONCEPT_ASSERT((RefEdgeListGraphConcept<CRefAdjList<TWeight> >));
    BOOST_CONCEPT_ASSERT((RefIncidenceGraphConcept<CRefAdjList<TWeight> >));
    BOOST_CONCEPT_ASSERT((RefAdjacencyGraphConcept<CRefAdjList<TWeight> >));
};


template class CRefEdgeListConceptCheck<no_edge_weight>;
template class CRefEdgeListConceptCheck<int>;
template class CRefEdgeListConceptCheck<double>;

template class CRefAdjListConceptCheck<no_edge_weight>;
template class CRefAdjListConceptCheck<int>;
template class CRefAdjListConceptCheck<double>;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
}

