function v = svm_dual_predict(S, K, alpha, b)
%SVM_DUAL_PREDICT SVM prediction based on dual formulation solution
%
%   v = SVM_DUAL_PREDICT(S, K, alpha);
%   v = SVM_DUAL_PREDICT(S, K, alpha, b);
%
%       Evaluates the prediction values based on the solution to the
%       dual QP problem, using the following formula
%
%       For classification:
%
%           v = w' * x + b = sum_i y_i alpha_i k(x_i, x) + b
%
%       Input arguments:
%       - S:        The SVM problem
%       - K:        The kernel matrix with respect to the samples on which
%                   the prediction is to be made
%
%                   Suppose there are N testing samples, then K will be
%                   a matrix of size S.n x N, where K(i, j) is the kernel
%                   value between the i-th training sample and the j-th
%                   testing sample.
%
%       - alpha:    The alpha coefficients
%       - b:        The offset value (if omitted, b is set to zero).
%
%       Output arguments:
%       - v:        The predicted values: size 1 x N.
%

% Created by Dahua Lin, on Jan 18, 2012


%% verify input

if ~is_svm_problem(S)
    error('svm_dual_predict:invalidarg', 'S should be a svm-problem struct.');
end

if ~(isfloat(K) && isreal(K) && ndims(K) == 2 && size(K,1) == S.n)
    error('svm_dual_predict:invalidarg', 'The kernel matrix K is invalid.');
end

if ~(isfloat(alpha) && isreal(alpha) && ndims(alpha) == 2 && size(alpha,2) == 1)
    error('svm_dual_predict:invalidarg', ...
        'alpha should be a real column vector.');
end

if nargin < 4
    b = 0;
else
    if ~(isfloat(b) && isreal(b) && isscalar(b))
        error('svm_dual_predict:invalidarg', 'b should be a real scalar.');
    end
end

%% main

switch S.type
    case 'class'
        v = (S.y .* alpha.') * K;
        
    otherwise
        error('svm_dual_predict:invalidarg', ...
            'Unsupported SVM problem type %s for svm_dual_predict', S.type);
end

if ~isequal(b, 0)
    v = v + b;
end

