/***************************************************************
 *
 *  gr_flood_cimp.cpp
 *
 *  The C++ mex implementation of gr_dijkstra_cimp
 *
 *  Created by Dahua Lin, on Jan 26, 2012
 *
 ***************************************************************/


#include "../../clib/smi_graph_mex.h"
#include <bcslib/array/array1d.h>
#include <bcslib/graph/graph_traversal.h>
#include <vector>

using namespace bcs;
using namespace bcs::matlab;
using namespace smi;

struct whole_roi
{
    bool viable(const vertex_t& v)
    {
        return true;
    }
};

struct mask_roi
{
    const bool *msk;
    
    mask_roi(const bool *msk_) : msk(msk_)
    {
    }
    
    bool viable(const vertex_t& v)
    {
        return msk[v.index()];
    }
};

template<class ROI>
struct bfs_vis : public trivial_traversal_agent<ginclist_t>
{
    std::vector<gint>& vs;
    ROI& roi;
    
    bfs_vis(std::vector<gint>& vs_, ROI& roi_) 
    : vs(vs_), roi(roi_) { }
    
    bool discover(const vertex_t& u, const vertex_t& v)
    {
        vs.push_back(v.id);
        return true;
    }
    
    bool examine(const vertex_t& u, const vertex_t& v, gvisit_status)
    {
        return roi.viable(v);
    }
};


template<class ROI>
struct bfs_visd : public trivial_traversal_agent<ginclist_t>
{            
    std::vector<gint>& vs;
    ROI& roi;
    array1d<int32_t>& dmap;
    
    bfs_visd(std::vector<gint>& vs_, ROI& roi_, array1d<int32_t>& dmap_) 
    : vs(vs_), roi(roi_), dmap(dmap_) { }
    
    bool discover(const vertex_t& u, const vertex_t& v)
    {
        vs.push_back(v.id);
        dmap[v.index()] = dmap[u.index()] + 1;
        return true;
    }
    
    bool examine(const vertex_t& u, const vertex_t& v, gvisit_status)
    {
        return roi.viable(v);
    } 
};



void do_flood(const ginclist_t& G, const vertex_t *s, gint ns, const bool *msk,
        int nlhs, mxArray *plhs[])
{
    if (nlhs <= 1)
    {
        std::vector<gint> vs;
        
        if (!msk)
        {
            whole_roi roi;
            bfs_vis<whole_roi> vis(vs, roi);
            breadth_first_traverse(G, vis, s, s + ns);        
        }
        else
        {
            mask_roi roi(msk);
            bfs_vis<mask_roi> vis(vs, roi);
            breadth_first_traverse(G, vis, s, s + ns); 
        }
        
        plhs[0] = to_matlab_row(vs).mx_ptr();
    }
    else
    {
        std::vector<gint> vs;
        array1d<int32_t> dmap(G.nvertices(), int32_t(0));
                
        if (!msk)
        {
            whole_roi roi;
            bfs_visd<whole_roi> vis(vs, roi, dmap);
            breadth_first_traverse(G, vis, s, s + ns);        
        }
        else
        {
            mask_roi roi(msk);
            bfs_visd<mask_roi> vis(vs, roi, dmap);
            breadth_first_traverse(G, vis, s, s + ns); 
        }
        
        std::vector<int32_t> ds;
        for (std::vector<int32_t>::const_iterator it = vs.begin();
            it != vs.end(); ++it)
        {
            ds.push_back(dmap[(*it) - 1]);
        }
        
        plhs[0] = to_matlab_row(vs).mx_ptr();
        plhs[1] = to_matlab_row(ds).mx_ptr();
    }
}



/**
 * Inputs:
 *
 *  [0] G:      The graph
 *  [1] s:      The vector of sources
 *  [2] msk:    The mask of ROI
 *
 * Outputs:
 *  [0] vs:     visited vertices
 *  [1] ds:     distances to sources 
 *
 */
void bcsmex_main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // take inputs
    
    const_marray mG(prhs[0]);
    ginclist_t G = mx_to_ginclist(mG);
    
    const_marray mS(prhs[1]);
    gint ns = mS.nelems(); 
    const vertex_t* s = mS.data<vertex_t>();
    
    const_marray mMsk(prhs[2]);
    const bool *msk = mMsk.is_empty() ? (const bool*)(0) : mMsk.data<bool>();
    
    // solve
    
    do_flood(G, s, ns, msk, nlhs, plhs);
    
}

BCSMEX_MAINDEF

