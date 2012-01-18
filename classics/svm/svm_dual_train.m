function [R, alpha] = svm_dual_train(S, K, solver)
%SVM_DUAL_TRAIN SVM training based on dual QP formulation
%
%   R = SVM_DUAL_TRAIN(S);
%   R = SVM_DUAL_TRAIN(S, K);
%   R = SVM_DUAL_TRAIN(S, [], solver);
%   R = SVM_DUAL_TRAIN(S, K, solver);
%
%       Trains a support vector machine (SVM) based on the dual QP
%       formulation of the problem given in S.
%
%       Input arguments:
%       - S:        The struct representing the SVM problem
%
%       - K:        The pre-computed kernel matrix. If not available,
%                   one can input K as an empty array.
%
%       - solver:   User-supplied solver (in form of a function handle).
%                   If omitted, the default solver that invokes MATLAB
%                   quadprog will be used.
%
%       In output, R is a struct that represents the solution, which 
%       has the following fields:
%       
%       - tag:      a fixed string: 'svm-dual-sol';
%       - type:     The problem type;
%       - n:        The number of support vectors
%       - X:        The matrix comprised of support vectors [d x n]
%       - a:        The coefficient on the support vectors. [1 x n]   
%       - b:        The offset scalar
%       - kernel_type
%       - kernel_params
%
%   [R, alpha] = SVM_DUAL_TRAIN( ... );
%
%       additionally returns the original QP solution alpha.
%
%   With the solution R, the prediction on a given set of samples X can 
%   be made by
%
%       v = R.a * svm_kernel(R, X) + R.b;
%
%   This can also be accomplished by invoking svm_dual_predict.
%

% Created by Dahua Lin, on Jan 18, 2012
%

%% process inputs

if nargin < 2
    K = [];
end

if nargin < 3
    solver = svm_default_solver();
else
    if ~isa(solver, 'function_handle')
        error('svm_dual_train:invalidarg', ...
            'solver must be a function handle.');
    end
end

%% main

% construct and solve QP

[P, K] = svm_dual(S, K);
alpha = solver(P);

% select support vectors

rtol = 1e-8;

switch S.type
    case 'class'        
        si = find(alpha > rtol * max(alpha));
        a = S.y(si) .* alpha(si).';
        
    case 'regress'
        n = S.n;
        a1 = alpha(1:n);
        a2 = alpha(n+1:2*n);
        a = a1 - a2;
        
        aa = abs(a);        
        si = find(aa > rtol * max(aa));
        a = a(si).';
        
    case 'rank'
        n = S.n;
        [I, J] = find(S.G);
        a1 = aggreg(alpha, n, I, 'sum');
        a2 = aggreg(alpha, n, J, 'sum');
        a = a1 - a2;
        
        aa = abs(a);
        si = find(aa > rtol * max(aa));
        a = a(si).';        
end

% solve offset (b)

b = svm_dual_offset(S, K, alpha, si);

% make solution

R.tag = 'svm-dual-sol';
R.type = S.type;
R.n = numel(si);
R.X = S.X(:, si);
R.a = a;
R.b = b;
R.kernel_type = S.kernel_type;
R.kernel_params = S.kernel_params;




