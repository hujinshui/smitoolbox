function [L, M, info] = kmeans_ex(X, M0, varargin)
% Extended K-means algorithm
%
%   [L, M] = kmeans_ex(X, K, ...);
%   [L, M] = kmeans_ex(X, M0, ...);
%       
%       This function implements the K-means algorithm, which is an
%       extension of the standard implementation.
%
%       Input: 
%       - X:    The matrix of input samples. Each column in X corresponds
%               to one sample.
%       - K:    The number of clusters. 
%       - M0:   The initial centers.
%   
%       Note that if K instead of M0 is input, the function will invoke
%       a method to initialize the initial centers. The option 'init' 
%       controls how to do the initialization. By default, it uses km++ 
%       (Kmeans++) method.
%
%       Output:
%       - L:    The assignment of samples to clusters. Suppose X contains
%               n samples, then L is a 1 x n row vector, whose values are
%               in {1, 2, ..., K}, and L(i) corresponds to X(:,i).
%       - M:    The resultant means, which is a matrix of size d x K,
%               where d is the space dimension, and K is the number of
%               clusters.
%
%       Options:
%       
%       The user can specify options (in name/value pairs or struct) to 
%       control the algorithm. The options are listed as below:
%
%       - max_iter: the maximum number of iterations (default = 100)
%
%       - tol_c:    the maximum allowable number of label changes at
%                   convergence. (default = 0)
%
%       - display:  the level of information displaying (default = 'off')
%                   It can take either of the following values:
%                   - 'off':    display nothing
%                   - 'iter':   display information for each iteration
%                   - 'final':  display final result information
%
%       - uc_warn:  whether to raise a warning if not converged.
%                   (default = false).
%
%       - init:     the method to initialize centers. It can take either
%                   of the following values (default = 'km++'):
%                   - 'random':  randomly pick K distinct samples as centers
%                   - 'km++':    randomly pick K samples using Kmean++
%                   - 'mdinit':  randomly select the first one, and then
%                                sequentially select furthest samples 
%                                from selected centers.                   
%
%       - on_nil:   the action taken to deal with "nil clusters", the 
%                   cluster to which no sample is assigned.
%                   (default = 'repick++').
%                   It can take either of the following values:
%                   - 'repick':     randomly repick a sample as new center
%                   - 'repick++':   randomly repick a sample as new center
%                                   using moderated scheme as in Kmeans++
%                   - 'mdpick':     pick the sample that is farthest to
%                                   all other centers
%                   - 'error':      raise an error.
%                   - 'keep':       simply keep that center without change.  
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%     
%
%   Remarks
%   -------
%       Briefly, compared to conventional implementation, this function
%       has additional features, as follows:
%
%       - It allows the use to choose between conventional implementation
%         and accelerated implementation. 
%
%         The conventional way computes the distances between all pairs of
%         samples and centers at every iteration, which leads to a great
%         deal of redundancy. The accelerated way reduces such redundant
%         computation exploiting triangle inequality, as described by the
%         following paper:
%
%           Charles Elkan. "Using the Triangle Inequality to Accelerate 
%           k-Means". Proceedings of 20th International Conference on
%           Machine Learning (ICML-03), Washington DC. 2003.
%
%         The basic idea of the accelerated method is to keep track of 
%         lower/upper bounds of relevant distances, and use triangle
%         inequality to identify unneccessary distance computation.
%         
%         We note that while the accelerated way avoids unnecessary 
%         computation, it would on the other hand incurs overhead of
%         book-keeping and subset selection. Generally, we recommend 
%         using the accelerated way when K is not very small.
%
%         Moreover, in our implementation, we avoid necessary
%         re-computation of mean-vectors, by only computing those affected
%         ones at each iteration.
%
%       - It allows the use of user-defined distance. The user can
%         implement its own distance by setting the dist_func option.
%         The function should support the following syntax:
%
%               R = dist_func(X, Y, 'd');
%                   computes the pairwise distances between X and Y.
%
%               R = dist_func(X, Y, 'c');
%                   computes the pairwise costs between X and Y. 
%                   Note that the cost used in objective function is not 
%                   necessary the same as distance. For example, in the
%                   standard Euclidean case, the cost of the square of
%                   the distance.
%
%               R = dist_func(X, [], 'm');
%               R = dist_func(X, w, 'm');
%                   computes the (weighted) mean of X.
%
%               R = dist_func(D, [], 't');
%                   computes the costs based on distances.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%
    
    






