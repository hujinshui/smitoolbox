function v = svm_dual_predict(R, X)
%SVM_DUAL_PREDICT SVM prediction based on dual formulation solution
%
%   v = SVM_DUAL_PREDICT(R, X);
%
%       Evaluates the prediction values based on the solution to the
%       dual QP problem, using the following formula
%
%           v = w' * x + b 
%             = sum_i a_i * k(sv_i, x) + b
%
%       Input arguments:
%       - R:        The svm-dual-sol struct that captures the SVM dual 
%                   problem solution.
%
%       - X:        The testing samples upon which the predicted values
%                   are to be evaluated.
%
%       Output arguments:
%       - v:        The predicted values: size 1 x N.
%

% Created by Dahua Lin, on Jan 18, 2012


%% main

v = R.a * svm_kernel(R, X);
b = R.b;

if ~isequal(b, 0)
    v = v + b;
end

