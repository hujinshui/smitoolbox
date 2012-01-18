function S = svm_problem(type, X, Y, C, kernel, tol)
%SVM_PROBLEM Support Vector Machine (SVM) Learning Problem
%
%   S = SVM_PROBLEM('class', X, y, C);
%   S = SVM_PROBLEM('class', X, y, C, kernel);
%
%       Constructs the problem struct for learning a binary SVM classifier.
%
%       The primal formulation is given by
%
%           minimize (1/2) * ||w||^2 + sum_i C_i * xi_i
%
%               with xi_i = max(0, 1 - y_i * (w' * x_i + b))
%
%       Here, if kernel is omitted, the problem is to learn a linear SVM,
%       otherwise, it is to learn a kernel SVM.
%
%
%   S = SVM_PROBLEM('regress', X, y, C, [], tol);
%   S = SVM_PROBLEM('regress', X, y, C, kernel, tol);
%
%       Constructs the problem struct for SVM regression.
%
%       The primal formulation is given by
%
%           minimize (1/2) * ||w||^2 + sum_i C_i * xi_i
%
%               with xi_i = max(0, abs(w' * x_i + b - y_i) - tol).
%
%
%   S = SVM_PROBLEM('rank', X, G);
%   S = SVM_PROBLEM('rank', X, G, [], kernel);
%
%       Constructs the problem struct for SVM ranking.
%
%       The primal formulation is given by
%
%           minimize (1/2) * ||w||^2 + sum_{ij} G_{ij} xi_{ij}
%
%               with xi_{ij} = max(0, 1 - w' * (x_i - x_j))
%
%   
%   Input arguments:
%       - X:        The feature matrix, of size d x n.
%
%       - y:        The vector of labels (for classification) or 
%                   responses (for ranking), of size 1 x n.
%                   For multi-class problem, the number of classes m
%                   are set to max(y).
%
%       - C:        The weights of slack variables, which can be 
%                   either a scalar, or a vector of size 1 x n to 
%                   assign sample-specific weights.
%
%       - kernel:   The kernel, which can be in either of the following
%                   form:
%
%                   - 'linear':         Linear kernel: k(x,y) = x' * y.
%                   - {'gauss', sigma}: Gaussian kernel: 
%                           k(x,y) = exp(-||x - y||^2/(2*sigma^2)).
%                   - {'poly', [a, b, d]}: Polynomial kernel:
%                           k(x,y) = (a * (x'*y) + b)^d
%                   - {'invmq', [a, b]}:  Inverse multiquadratic kernel:
%                           k(x,y) = 1 / (a*||x-y||^2 + b)
%                   - {'sigmoid', [a, b]}: Sigmoid kernel:
%                           k(x,y) = tanh(a * (x'*y) + b)
%
%                   Or it can be a function handle, such that
%                   kernel(X, Y) returns an m x n matrix K with
%                   K(i, j) being the kernel value of X(:,i) and Y(:,j).
%
%                   If kernel is omitted, then the linear kernel is
%                   assumed.
%       
%       Note that if a pre-computed kernel is provided, the sample matrix
%       X can be input as an empty array.
%
%   The output S is a struct with the following fields:
%
%       - tag:      A fixed string: 'svm-problem'.
%
%       - type:     The problem type, whose value can be either of:
%                   'class', 'multi-class', 'regress', and 'rank'.
%
%       - n:        The number of training samples
%
%       - m:        The number of classes. 
%
%       - X:        The sample matrix.
%
%       - y:        The vector of labels or responses.
%
%       - G:        The ranking coefficient matrix. (G is empty, if 
%                   S is not a ranking problem)
%
%       - C:        The slack weight (a scalar or a 1 x n vector).
%
%       - tol:      The tolerance for regression
%
%       - kernel_type:  The type of kernel, which can be either of:
%                       'linear', 'gauss', 'poly', 'invmq', or 'custom'.
%                       ('custom' means using user-supplied function).
%
%       - kernel_params:  depends on the kernel type:
%                         - for linear kernel, it is empty.
%                         - for gauss kernel, it has a field sigma
%                         - for poly kernel, it has fields: a, b, d
%                         - for invmq kernel, it has fields: a, b
%                         - for sigmoid kernel, it has fields: a, b
%                         - for custom kernel, it is the function handle.                           
%

% Created by Dahua Lin, on Jan 17, 2012
%

%% verify inputs

if ~ischar(type)
    error('svm_problem:invalidarg', 'The SVM type should be a string.');
end

if ~(isnumeric(X) && ndims(X) == 2)
    error('svm_problem:invalidarg', 'X should be a numeric matrix.');
end
n = size(X, 2);

switch type
    case 'class'
        if isfloat(Y) && isreal(Y) && isequal(size(Y), [1 n])
            m = 2;
            y = Y;
            G = [];
        else
            error('svm_problem:invalidarg', 'y should be an 1 x n real vector.');
        end
        
    case 'regress'                
        if isfloat(Y) && isreal(Y) && isequal(size(Y), [1 n])
            m = 1;
            y = Y;
            G = [];
        else
            error('svm_problem:invalidarg', 'y should be an 1 x n real vector.');
        end
                
    case 'rank'
        if isfloat(Y) && isreal(Y) && isequal(size(Y), [n n])
            m = 1;
            y = [];
            G = Y;
        else
            error('svm_problem:invalidarg', 'G should be an n x n real matrix.');
        end
        
    otherwise
        error('svm_problem:invalidarg', 'The SVM type is invalid.');
end

if ~(isfloat(C) && ((isscalar(C) && C >= 0) || ...
        (isvector(C) && numel(C) == n)))
    error('svm_problem:invalidarg', ...
        'C should be either a positive scalar or a real vector of length n.');
end

if nargin < 5 || isempty(kernel)
    kertype = 'linear';
    kerparam = [];
else
    [kertype, kerparam] = check_kernel(kernel);  
end

if strcmp(type, 'regress')
    if ~(isfloat(tol) && isreal(tol) && isscalar(tol) && tol > 0)
        error('svm_problem:invalidarg', 'tol is invalid.');        
    end
else
    tol = [];
end

 
%% make output

S.tag = 'svm-problem';
S.type = type;
S.n = n;
S.m = m;
S.X = X;
S.y = y;
S.G = G;
S.C = C;
S.tol = tol;

S.kernel_type = kertype;
S.kernel_params = kerparam;

if strcmp(kertype, 'pre')
    S.kermat = kernel;
end


%% kernel checking

function [kt, kp] = check_kernel(kernel)
% verify kernels

if strcmpi(kernel, 'linear')
    kt = 'linear';
    kp = [];
elseif iscell(kernel)
    kt = lower(kernel{1});
    r = kernel{2};
    
    switch kt
        case 'gauss'
            if isfloat(r) && isscalar(r) && r > 0
                kp.sigma = r;
            else
                error('svm_problem:invalidarg', ...
                    'The parameter for the gauss kernel is invalid.');
            end
            
        case 'poly'                        
            if isfloat(r) && numel(r) == 3 && all(r >= 0)
                kp.a = r(1);
                kp.b = r(2);
                kp.d = r(3);
            else
                error('svm_problem:invalidarg', ...
                    'The parameter for the poly kernel is invalid.');
            end
            
        case 'invmq'                        
            if isfloat(r) && numel(r) == 2 && all(r >= 0)
                kp.a = r(1);
                kp.b = r(2);
            else
                error('svm_problem:invalidarg', ...
                    'The parameter for the invmq kernel is invalid.');
            end
            
        case 'sigmoid'                        
            if isfloat(r) && numel(r) == 2 && all(r >= 0)
                kp.a = r(1);
                kp.b = r(2);
            else
                error('svm_problem:invalidarg', ...
                    'The parameter for the sigmoid kernel is invalid.');
            end
            
    end
elseif isa(kernel, 'function_handle')
    kt = 'custom';
    kp = kernel;
else
    error('svm_problem:invalidarg', ...
        'The input kernel is invalid.');
end

