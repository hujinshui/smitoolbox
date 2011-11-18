function model = primal_svm_train(X, y, c, varargin)
% Train a linear or kernel SVM using primal training
%
%   model = primal_svm_train(X, y, c, ...);
%
%       This function learns an SVM from data by directly optimizing
%       the primal objective as an unconstrained problem, as follows.
%       
%       For linear SVM, the problem is formalized as follows:
%
%           minimize (1/2) ||w||^2 
%                    + sum_i c_i * L( y_i * f(w, b, x_i) )
%
%           with f(w, b, x) = w * x + b.
%
%       For kernel SVM with kernel function f, the problem is
%
%           minimize (1/2) beta' * K * beta 
%                    + sum_i c_i * L( y_i * f(beta, b, x_i) )
%
%           with f(beta, b, x) = sum_i beta_i * k(x_i, x) + b.
%
%       Input arguments:
%       - X:        the input feature matrix, with each column being 
%                   an sample. Suppose there are n samples of dimension d,
%                   X should be a matrix of d x n, with X(:,i)
%                   corresponding to the i-th sample.
%
%       - y:        the label vector, where y(i) is the label of X(:,i),
%                   whose value can be -1 or 1.
%
%       - c:        the weights of sample loss. It can be either a 
%                   scalar or a vector of length n.
%
%       Output arguments:
%       - model:    the trained model. Depending on whether kernel is 
%                   used, it returns an object of class kernel_svm or
%                   linear_svm.
%                   (To use kernel, one has to explicitly specify
%                    the kernel function in additional options).
%
%       Additional options (which can be given via name/value pairs)
%
%       - kernel:   the kernel function handle. It will be called to
%                   compute kernel when necessary in the following way:
%
%                   kernel(X, Y) should return a m x n matrix K when
%                   X and Y respectively have m and n columns, and 
%                   K(i, j) is the kernel value between X(:,i) and Y(:,j).
%
%                   Default = [], which indicates that linear SVM is to be
%                   trained.
%
%       - kermat:   the pre-computed kernel matrix of size n x n.
%                   Note that even when kermat is given, the kernel 
%                   function handle still need to be specified for
%                   constructing a generalizable classifier.
%
%                   Default = [], which indicates that kernel matrix is
%                   not pre-computed.
%
%       - schedule: For sample schedule for multiple-stage training. 
%                   In particular, schedule should be a cell array of 
%                   index vectors. If this is specified, in the first 
%                   stage, schedule{1} is used to train the model at
%                   first, and then schedule{2} is used in the second
%                   stage, with the output of schedule{1} as initial guess
%                   of the solution.
%
%                   One can set schedule{k} to [] to indicate using 
%                   the whole sample set in the k-th stage.
%
%                   Default = {[]}, which indicates that there is only one
%                   stage that takes in all samples as input.
%
%       - initsol:  The initial guess of the solution, which is a 
%                   vertical concatenation as [w;b].
%
%       - huberwidth:   The width of transition band in huber function.
%                       Default = 1e-3.
%
%       - regb:     The regularization coefficient of b. Default = 1e-8.
%
%       - verbose:  whether to show the training progress. Default = true.
%
%       - optimopt: The option for the optimization routine (newtonfmin).
%                   It can be obtained using smi_optimset. 
%                   Default = newtonfmin('options');
%       

%   History
%   -------
%       - Created by Dahua Lin, on April 21, 2011
%


%% verify input argument

% data

if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
    error('primal_svm_train:invalidarg', 'X should be a real matrix.');
end
if ~isa(X, 'double'); X = double(X); end

[d, n] = size(X);

if ~(isnumeric(y) && isvector(y) && length(y) == n && isreal(y))
    error('primal_svm_train:invalidarg', 'y should be a real vector of length n.');    
end
if ~isa(y, 'double'); y = double(y); end
if size(y, 1) > 1; y = y.'; end

if ~(isfloat(c) && isreal(c) && (isscalar(c) || (isvector(c) && length(c) == n)))
    error('primal_svm_train:invalidarg', ...
        'c should be a scalar or a real vector of length n.');
end
if ~isa(c, 'double'); c = double(c); end
if size(c, 1) > 1; c = c.'; end


% options

kernel = [];
kermat = [];
schedule = {[]};
initsol = [];
huberwidth = 1e-3;
regb = 1e-8;
optimopt = [];
verbose = false;


if ~isempty(varargin)
    
    onames = varargin(1:2:end);
    ovals = varargin(2:2:end);
    
    if ~(length(onames) == length(ovals) && iscellstr(onames))
        error('primal_svm_train:invalidarg', ...
            'The options are not correctly specified.');
    end
    
    for i = 1 : numel(onames)        
        cn = onames{i};
        cv = ovals{i};
        
        switch lower(cn)            
            case 'kernel'
                if ~isa(cv, 'function_handle')
                    error('primal_svm_train:invalidarg', ...
                        'The kernel should be a function handle.');
                end
                kernel = cv;
                    
            case 'kermat'
                if ~(isfloat(cv) && isequal(size(cv), [n n]))
                    error('primal_svm_train:invalidarg', ...
                        'The kermat should be a numeric matrix of size n x n.');
                end
                kermat = cv;
                
            case 'schedule'
                if ~(iscell(cv) && all(cellfun(@isnumeric, cv(:))))
                    error('primal_svm_train:invalidarg', ...
                        'The schedule should be a cell array of index arrays.');
                end
                schedule = cv;
                
            case 'initsol'
                if ~(isfloat(cv) && size(cv, 2) == 1 && ndims(cv) == 2 && isreal(cv))
                    error('primal_svm_train:invalidarg', ...
                        'The initsol should be a real column vector.');
                end
                initsol = cv;                
                
            case 'huberwidth'
                if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv > 0 && cv < 1)
                    error('primal_svm_train:invalidarg', ...
                        'huberwidth should be a real value in (0, 1).');
                end
                huberwidth = cv;
                
            case 'regb'
                if ~(isfloat(cv) && isscalar(cv) && isreal(cv) && cv >= 0)
                    error('primal_svm_train:invalidarg', ...
                        'regb should be a non-negative real scalar.');
                end
                regb = cv;
                
            case 'optimopt'
                if ~(isstruct(cv) && isscalar(cv))
                    error('primal_svm_train:invalidarg', ...
                        'optimopt should be a struct.');
                end
                optimopt = cv;
                
            case 'verbose'
                if ~((isnumeric(cv) || islogical(cv)) && isscalar(cv))
                    error('primal_svm_train:invalidarg', ...
                        'verbose should be a numeric/logical scalar.');
                end
                verbose = logical(cv);                 
                
            otherwise
                error('primal_svm_train:invalidarg', ...
                    'Unknown option name %s', cn);
        end                
    end    
end

use_kernel = ~isempty(kernel);

if ~isempty(kermat)
    if ~use_kernel
        error('primal_svm_train:invalidarg', ...
            'kernel must be specified when kermat is given.');
    end
end

if use_kernel
    dsol = n + 1;
else
    dsol = d + 1;
end
if isempty(initsol)
    initsol = zeros(dsol, 1);
else
    if size(initsol, 1) ~= dsol
        error('primal_svm_train:invalidarg', ...
            'The dimension of initsol is incorrect, which should be %d.', dsol);
    end
end

if isempty(optimopt)
    optimopt = newtonfmin('options');
end
optimopt.DirectNewton = true;


%% main

if verbose
    if use_kernel
        fprintf('Training Primal kernel-SVM on %d samples ...\n', n);
    else
        fprintf('Training Primal linear-SVM on %d samples of dim %d ...\n', n, d);
    end
end

% prepare kernel matrix

if use_kernel
    if verbose
        fprintf('Preparing kernel matrix ...\n');
    end
    
    if isempty(kermat)
        K = kernel(X, X);
    else
        K = kermat;
    end
end

% initialize

sch = [];
sol = initsol;
lof = @(z) svm_huber_loss(z, huberwidth);


% main loop

nsch = numel(schedule);
for t = 1 : nsch
    
    pre_sch = sch;
    sch = schedule{t};
    
    if verbose
        if isempty(sch)
            fprintf('Training Stage %d: on %d samples ...\n', t, n);
        else
            fprintf('Training Stage %d: on %d samples ...\n', t, numel(sch));
        end
    end    
    
    if use_kernel
        
        % expand sol
        if ~isequal(sch, pre_sch)
            if isempty(pre_sch)
                temp_sol = sol;
            else
                temp_sol = zeros(n+1, 1);
                temp_sol(pre_sch) = sol(1:end-1);
                temp_sol(n+1) = sol(end);
            end
            sol = [temp_sol(sch); temp_sol(n+1)];
        end        
        
        % make Q and F
        if isempty(sch)
            Q = K;
            yc = y;
        else
            Q = K(sch, sch);
            yc = y(sch);
        end
        F = [];
        
    else % not use kernel
        
        % make Q and F
        Q = 1;
        if isempty(sch)
            F = X;
            yc = y;
        else
            F = X(:, sch);
            yc = y(sch);
        end
    end
        
    fd = @(s) primal_svm_obj(s, Q, F, yc, c, lof, regb, 'd');
    
    sol = newtonfmin(fd, sol, optimopt);
    
end


% make output

if use_kernel
    beta = sol(1:n);
    b = sol(n+1);    
    model = kernel_svm.from_psol(X, y, beta, b, kernel);
else
    model = linear_svm.from_sol(X, y, sol);
end







