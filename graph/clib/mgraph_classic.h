/********************************************************************
 *  
 *  mgraph_classic.h
 *
 *  The header files for some classic graph theory algorithms
 *
 *  Created by Dahua Lin, on Oct 1, 2009
 *
 ********************************************************************/

#ifndef SMI_MGRAPH_CLASSIC_H
#define SMI_MGRAPH_CLASSIC_H

#include "mgraph.h"

namespace smi
{
 
/**
 * Performs breadth-first traverse from a specified node
 *
 * @param G the (outgoing) neighborhood of the graph
 * @param v0 the index of the starting node (inclusive)
 * @param r the pre-allocated buffer to store traversed nodes
 * @param visited the bool array indicating which nodes have been visited 
 *        (should be properly initialized/set before)
 *
 * @return the number of traversed nodes
 */   
int breadth_first_traverse(const GNeighborHood& G, int v0, int *r, bool *visited);

/**
 * Performs depth-first traverse from a specified node
 *
 * @param G the (outgoing) neighborhood of the graph
 * @param v0 the index of the starting node (inclusive)
 * @param r the pre-allocated buffer to store traversed nodes 
 * @param visited the bool array indicating which nodes have been visited 
 *        (should be properly initialized/set before)
 *
 * @return the number of traversed nodes
 */
int depth_first_traverse(const GNeighborHood& G, int v0, int *r, bool *visited);
    
}


#endif

