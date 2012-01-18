function [alpha, b] = svm_dual_train(S, K, solver)
%SVM_DUAL_TRAIN SVM training based on dual QP formulation
%
%   [alpha, b] = SVM_DUAL_TRAIN(S);
%   [alpha, b] = SVM_DUAL_TRAIN(S, K);
%   [alpha, b] = SVM_DUAL_TRAIN(S, [], solver);
%   [alpha, b] = SVM_DUAL_TRAIN(S, K, solver);
%
%       Trains a support vector machine (SVM) based on the dual QP
%       formulation of the problem given in S.
%
%       The solution comprises the alpha coefficients and the offset
%       scalar b.
%
%       The caller can supply a pre-computed kernel matrix K to 
%       facilitate the computation if available. 
%
%       The caller can also use a customized solver by providing it in 
%       form of a function handle. Otherwise, the default solver will
%       be used.
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

[P, K] = svm_dual(S, K);
alpha = solver(P);
b = svm_dual_offset(S, K, alpha);

