function P = kernel_svm_prob(K, y, c)
% Construct the dual QP problem for standard kernel SVM
%
%   The standard kernel support vector machine problem (in dual form) 
%   is formulated as follows
%
%       minimize (1/2) * a' * H * a - sum(a)
%
%       s.t. 0 <= a_i <= c
%            sum_i a_i y_i = 0.
%
%       Here, H is an n x n matrix with 
%
%       H(i,j) = y_i * y_j * K(i, j)
%
%       and K(i, j) is the kernel value of x_i and x_j.
%
%   
%   P = kernel_svm_prob(K, y, c);
%
%       constructs a QP problem P as formulated above and returns it.
%
%       Input arguments:
%   
%       - K:    the kernel matrix of the input samples. Suppose there are
%               n samples, then the size of K is n x n. K can be either
%               a full matrix or a sparse matrix.
%
%       - y:    the label vector of length n. In particular, y(i) is
%               the label for X(i), whose value can be either -1 or 1.
%
%       - c:    the upper bound of a_i values.
%
%       The output P is a struct representing a QP problem. One can use
%       any available QP problem solver to solve it, and thus obtain
%       the alpha vector.
%       
%               

%   History
%   -------
%       - Created by Dahua Lin, on Apr 7, 2011
%

%% verify input arguments

n = size(K, 1);
if ~(isfloat(K) && ndims(K) == 2 && isreal(K) && ~isempty(K) && size(K,2) == n)
    error('kernel_svm_prob:invalidarg', ...
        'K should be a non-empty real square matrix.');
end
if ~isa(K, 'double'); K = double(K); end
    
     
if ~(isnumeric(y) && isvector(y) && isreal(y) && length(y) == n)
    error('kernel_svm_prob:invalidarg', ...
        'y should be a real vector of length n.');
end
if ~isa(y, 'double'); y = double(y); end
if size(y, 1) > 1; y = y.'; end   % turn into a row vector

if ~(isfloat(c) && isreal(c) && isscalar(c) && c > 0)
    error('kernel_svm_prob:invalidarg', ...
        'c should be a positive real scalar.');
end

%% main

H = bsxfun(@times, bsxfun(@times, K, y), y');
f = -ones(n, 1);

P = qp_problem(H, f, [], [], y, 0, 0, c);



