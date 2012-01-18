function K = svm_kernel(S, X)
%SVM_KERNEL SVM Kernel evaluation
%
%   K = SVM_KERNEL(S);
%
%       Evaluates the kernel matrix (i.e. Gram matrix) for the given
%       SVM problem or solution.
%
%       Here, S should be either a svm-problem or a svm-dual-sol struct.
%
%       In the output, K is a matrix of size S.n x S.n.
%
%   K = SVM_KERNEL(S, X);
%
%       Evaluates the kernel matrix between S.X and X. The output 
%       matrix K is of size S.n x size(X,2).
%

% Created by Dahua Lin, on Jan 18, 2011
%

%% verify input arguments

if ~(is_svm_problem(S) || is_svm_dual_sol(S))
    error('svm_kernel:invalidarg', ...
        'S should be a svm-problem or svm-dual-sol struct.');
end

if nargin < 2
    X = [];
else
    d = size(S.X, 1);
    if ~(isnumeric(X) && ndims(X) == 2 && size(X,1) == d)
        error('svm_kernel:invalidarg', 'The input matrix X is invalid.');
    end
end

%% main

kp = S.kernel_params;

switch S.kernel_type
    case 'linear'
        if isempty(X)
            X = S.X;
            K = X' * X;
        else
            K = S.X' * X;
        end
        
    case 'gauss'
        K = gausskernel(S.X, X, kp.sigma);
        
    case 'poly'
        K = polykernel(S.X, X, kp.a, kp.b, kp.d);
        
    case 'invmq'
        K = invmqkernel(S.X, X, kp.a, kp.b);
        
    case 'sigmoid'
        K = sigmoidkernel(S.X, X, kp.a, kp.b);
        
    otherwise
        error('svm_kernel:invalidarg', ...
            'Unsupported kernel type %s', S.kernel_type);
end
        



