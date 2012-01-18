function [w, b, xi] = svm_primal_train(S, solver)
%SVM_PRIMAL_TRAIN Train Support Vector Machine using Primal formulation
%
%   [w, b] = SVM_PRIMAL_TRAIN(S);
%   [w, b] = SVM_PRIMAL_TRAIN(S, solver);
%   [w, b, xi] = SVM_PRIMAL_TRAIN( ... );
%
%       trains a linear SVM using primal QP formulation.
%
%       Input arguments:
%       - S:        The SVM problem struct. 
%                   Note S.kernel_type must be 'linear'.
%
%       - solver:   The solver function handle. 
%                   If omitted, the function by default uses mstd_solve,
%                   which invokes MATLAB's built-in interior-point solver.
%
%       Output arguments:
%       - w:        The weight vector [d x 1].
%       - b:        The offset scalar.
%       - xi:       The vector of slack variables
%
%       The linear prediction is then given by w'*x+b.
%

% Created by Dahua Lin, on Jan 18, 2012
%

%% Prepare solver

if nargin < 2
    solver = svm_default_solver();
else
    if ~isa(solver, 'function_handle')
        error('svm_primal_train:invalidarg', ...
            'solver must be a function handle.');
    end
end

%% main

P = svm_primal(S);
sol = solver(P);

d = size(S.X, 1);

switch S.type
    case {'class', 'regress'}
        w = sol(1:d);
        b = sol(d+1);
        if nargout >= 3
            xi = sol(d+2:end);
        end
        
    case 'rank'
        w = sol(1:d);
        b = 0;
        if nargout >= 3
            xi = sol(d+1:end);
        end
end


