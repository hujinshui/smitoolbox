function [mu, Sw, Sb] = scattermat(X, K, L)
% Compute scatter matrices for labeled Data
%   
%   [mu, Sw, Sb] = scattermat(X, K, L);
%   [mu, Sw, Sb] = scattermat(X, K, L, w);
%
%       computes within-class (and between-class) scatter matrices
%       from labeled data.
%
%       Inputs:
%       - X:    a d x n matrix, with each column being sample.
%       - K:    the number of classes
%       - L:    the 1 x n label vector. L(i) is the label of sample X(:,i),
%               whose value should be in {1, 2, ..., K}.
%       - w:    the weights of samples (a vector of size 1 x n). 
%               If omitted, all samples have the same weight.
%
%       Outputs:
%       - mu:   the means of the samples
%       - Sw:   the within-class scatter matrix (i.e. pooled covariance)
%       - Sb:   the between-class scatter matrix (i.e. covariance of mean)
%

% Created by Dahua Lin, on Nov 22, 2010.
%

