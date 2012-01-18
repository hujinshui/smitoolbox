function [P, K] = svm_dual(S, K)
%SVM_DUAL Dual QP Problem of SVM Learning
%
%   P = svm_dual(S);
%       
%       constructs a QP problem corresponding to the dual formulation
%       of the SVM problem given by S.
%
%   P = svm_dual(S, K);
%
%       performs the construction using the pre-computed kernel matrix
%       K as the kernel matrix. If K is not provided, it will be computed
%       by the function.
%
%   [P, K] = svm_dual( ... );
%
%       additionally returns the kernel matrix.
%
%   Note that the solution to a dual problem is the vector of the alpha
%   coefficients.
%

% Created by Dahua Lin, on Jan 18, 2012
%

%% verify inputs

if ~is_svm_problem(S)
    error('svm_dual:invalidarg', 'S should be a svm-problem struct.');
end
n = S.n;

if nargin < 2 || isempty(K)
    K = [];
else
    if ~(isfloat(K) && isreal(K) && isequal(size(K), [n n]))
        error('svm_dual:invalidarg', 'K should be an n x n real matrix.');
    end
end

%% main

% get kernel matrix

if isempty(K)
    K = svm_kernel(S);
end

% construct problem

switch S.type
    case 'class'
        sdim = n;
        y = S.y;
        H = (y' * y) .* K;
        f = -ones(n, 1);
        Aeq = y;
        beq = 0;
        c = S.C.';        
        
    otherwise
        error('svm_dual:invalidarg', ...
            'Unsupported SVM problem type %s', S.type);
end


if isscalar(c)
    c = c * ones(sdim, 1);
end
P = qp_problem(H, f, [], [], Aeq, beq, zeros(sdim, 1), c);



