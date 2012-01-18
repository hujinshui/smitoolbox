function S = svm_problem(type, X, Y, C, kernel)
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
%   S = SVM_PROBLEM('multi-class', X, y, C);
%   S = SVM_PROBLEM('multi-class', X, y, C, kernel);
%
%       Constructs the problem struct for learning a multi-class SVM
%       classifier.
%
%       The primal formulation is given by
%
%           minimize (1/2) * sum_{k=1}^m ||w_k||^2 +
%                    sum_i sum_{k != y_i} C_i * xi_{i,k}
%
%           with xi_{i,k} = max(0, 1 - (w_{y_i}' * x_i - w_k * x_i))
%
%
%   S = SVM_PROBLEM('regress', X, y, C);
%   S = SVM_PROBLEM('multi-class', X, y, C, kernel);
%
%       Constructs the problem struct for SVM regression.
%
%       The primal formulation is given by
%
%           minimize (1/2) * ||w||^2 + sum_i C_i * xi_i
%
%               with xi_i = max(0, abs(w' * x_i + b - y_i) - tol).
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
%               with xi_{ij} = max(0, 1 - w' * (y_i - y_j))
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
%
%                   Or it can be an n x n pre-computed kernel matrix.
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
%       - kernel_type:  The type of kernel, which can be either of:
%                       'linear', 'gauss', 'poly', 'invmq', or 'pre'.
%                       ('pre' means pre-computed).
%
%       - kermat:   The pre-computed kernel matrix (it is empty if not
%                   available).
%
%       - kernel_params:  Empty when the kernel is linear or pre-computed.
%                         It is a struct with fields sigma for gauss
%                         kernel. It is a struct with fields a, b, and d
%                         for poly kernel, and with fields a and b for 
%                         invmq kernel.
%
%

% Created by Dahua Lin, on Jan 17, 2012
%

%% verify inputs

if ~ischar(type)
    error('svm_problem:invalidarg', 'The SVM type should be a string.');
end

if ~isempty(X)
    if ~(isnumeric(X) && ndims(X) == 2)
        error('svm_problem:invalidarg', 'X should be a numeric matrix.');
    end
end

if nargin < 5
    kertype = 'linear';
    kerparam = [];
else
    if ~isnumeric(kernel)
        [kertype, kerparam] = check_kernel(kernel);        
    else
        n = size(kernel, 1);
        if ~(isfloat(kernel) && ndims(kernel) == 2 && size(kernel,2) == n)
            error('svm_problem:invalidarg', ...
                'Pre-computed kernel must be a real square matrix.');
        end
        if ~isempty(X)
            if size(X, 2) ~= n
                error('svm_problem:invalidarg', ...
                    'The size of X is inconsistent with the pre-computed kernel.');                
            end
        end
        
        kertype = 'pre';
        kerparam = [];
    end
end

if ~strcmp(kertype, 'pre')
    if isempty(X)
        error('svm_problem:invalidarg', ...
            'X must be non-empty, when pre-computed kernel is not provided.');
    end
    n = size(X, 2);
end
    
switch type
    case 'class'
        if isfloat(Y) && isreal(Y) && isequal(size(Y), [1 n])
            m = 2;
            y = Y;
            G = [];
        else
            error('svm_problem:invalidarg', 'y is invalid.');
        end
        
    case 'multi-class'
        if isnumeric(Y) && isreal(Y) && isequal(size(Y), [1 n])
            m = double(max(Y));
            y = Y;
            G = [];
        else
            error('svm_problem:invalidarg', 'y is invalid.');
        end
        
    case 'regress'                
        if isfloat(Y) && isreal(Y) && isequal(size(Y), [1 n])
            m = 1;
            y = Y;
            G = [];
        else
            error('svm_problem:invalidarg', 'y is invalid.');
        end
        
    case 'rank'
        if isfloat(Y) && isreal(Y) && isequal(size(Y), [n n])
            m = 1;
            y = [];
            G = Y;
        else
            error('svm_problem:invalidarg', 'G is invalid.');
        end
        
    otherwise
        error('svm_problem:invalidarg', 'The SVM type is invalid.');
end

if ~(isfloat(C) && ((isscalar(C) && C >= 0) || ...
        (isvector(C) && numel(C) == n)))
    error('svm_problem:invalidarg', ...
        'C should be either a positive scalar or a real vector of length n.');
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
                    'The parameter for the poly kernel is invalid.');
            end
    end
else
    error('svm_problem:invalidarg', ...
        'The input kernel is invalid.');
end

