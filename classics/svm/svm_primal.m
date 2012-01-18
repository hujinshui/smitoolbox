function P = svm_primal(S)
%SVM_PRIMAL Primal QP Problem of SVM Learning
%
%   P = SVM_PRIMAL(S);
%
%       constructs a QP problem corresponding to the primal formulation
%       of the SVM problem represented by S.
%
%       The solution composition for different types of problems are
%       different.
%
%       - 'class':          sol = [w; b; xi]
%       - 'regress':        sol = [w; b; xi]
%       - 'rank':           sol = [w; xi]
%
%       Note this function only applies to the problem using linear kernel.
%

% Created by Dahua Lin, on Jan 17, 2012
%

%% verify input

if ~is_svm_problem(S)
    error('svm_primal:invalidarg', 'S should be a svm-problem struct.');
end

if ~strcmp(S.kernel_type, 'linear')
    error('svm_primal:invalidarg', ...
        'S.kernel_type must be linear for svm_primal.');
end


%% main

n = S.n;
d = size(S.X, 1);

% generate constraints 

switch S.type
    case 'class'
        X = S.X;
        y = S.y;        
        A = [bsxfun(@times, y, X).', y.']; 
        b = ones(n, 1);
        c = S.C;       
        
    case 'regress'
        X = S.X.';
        y = S.y.';
        tol = S.tol;
        A = [-X, -ones(n, 1); X, ones(n, 1)];
        b = [-(y + tol); y - tol]; 
        c = S.C;
        if ~isscalar(c)
            c = [c c];
        end
        
    case 'rank'
        [I, J, c] = find(S.G);        
        X = S.X;
        A = (X(:, I) - X(:, J)).';
        
    otherwise
        error('svm_primal:invalidarg', ...
            'Unsupported SVM model type %s', S.type);
end

% generate constraint matrix

cdim = size(A, 1);
if issparse(A) || 5 * d < cdim
    A = [A, speye(cdim)];
else
    A = [A, eye(cdim)];
end
sdim = size(A, 2);

% quadratic coefficient matrix

wdim = d;
H = sparse(1:wdim, 1:wdim, 1, sdim, sdim);

% linear coefficient vector

f = zeros(sdim, 1);
f(end-cdim+1:end) = c;

% lower-bound

lb = zeros(sdim, 1);
lb(1:end-cdim) = -inf;

% make problem struct

P = qp_problem(H, f, -A, -b, [], [], lb, []);


