function W = lda_ex(X, K, L, varargin)
% Extended Multi-class Linear Discriminant Analysis
%
%   T = lda_ex(X, K, L, ...);
%       performs multi-class linear discriminant analysis. 
%
%       Suppose the training set comprises n samples on a d-dimensional
%       vector space. Then X should be a matrix of size d x n. Each
%       column of X is a sample.
%
%       K is the number of classes. L is a 1 x n row vector of the class 
%       labels. In particular, L(i) is the label for X(:,i), whose value 
%       should be in {1, 2, ..., K}, where K is the number of classes.
%
%       In the output, T is a transform matrix of size p x d (p < d), 
%       which can transform d-dimensional input vectors into 
%       p-dimensional discriminant features.
%
%       One can further specify some of the following options to 
%       control the algorithm in form of name/value pairs. 
%       For options that are not explicitly specified, the default
%       values will be used.
%
%       
