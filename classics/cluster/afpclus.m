function [M, L] = afpclus(S, varargin)
%AFPCLUS Clustering by Affinity Propagation
%
%   [M, L] = AFPCLUS(S, ...);
%
%       Performs Affinity propagation to cluster data samples.
%
%       Input arguments:
%       - S:            The similarity matrix, of size n x n.
%                       Here, n is the number of samples.
%
%                       Note: S(i,i) is the preference to make the 
%                       i-th sample to be one of the centers.
%
%
%       Output arguments:
%       - M:            The vector of center indices [1 x K]
%
%       - L:            The assignment vector [1 x n]. In particular, 
%                       L(i) = k indicates that the i-th sample is
%                       associated with the examplar M(k).
%                       
%       One can specify other options to control the procedure, in form
%       of name/value pairs:
%
%       - 'MaxIter':    The maximum number of iterations (default = 500)
%
%       - 'LearnRate':  The learning rate (default = 0.5)
%                       updated = (1 - rate) * old + rate * new.
%
%       - 'ExamIntv':   The interval of examination (default = 10)
%
%                       The assignment is examined per ExamIntv iterations.
%                       If the assignment within this period does not
%                       change, the function considers the procedure
%                       as converged.
%
%       - 'Display':    The level of displaying (default = 'off')
%                       - 'off':    display nothing
%                       - 'final':  display when the procedure finishes
%                       - 'phase':  display at each phase
%                                   (the examination marks the end of
%                                    a phase)
%

% Created by Dahua Lin, on Mar 28, 2012


%% parse input arguments

n = size(S, 1);
if ~(isfloat(S) && isreal(S) && ~issparse(S) && size(S, 2) == n)
    error('afpclus:invalidarg', 'S should be a real square matrix,');
end

[maxiter, lrate, exam_intv, displevel] = check_opts(varargin);


%% main

% initialize

A = zeros(n, n);
R = zeros(n, n);

% main loop

t = 0;
converged = false;

while ~converged && t <= maxiter
            
    % compute responsibility
    
    R_pre = R;    
    R = compute_r(S, A);        
    if lrate < 1
        R = (1 - lrate) * R_pre + lrate * R;
    end
    
    % compute availability
    
    A_pre = A;
    A = compute_a(R);        
    if lrate < 1
        A = (1 - lrate) * A_pre + lrate * A;  
    end
    
    % examine
    
    if t == 0
        M = get_centers(A, R);
        last_exam = 0;
        
    elseif t == last_exam + exam_intv
        M_last = M;
        M = get_centers(A, R);
        converged = isequal(M, M_last);
        last_exam = t;
        
        if displevel >= 2
            ch = nnz(get_ass(S, M) - get_ass(S, M_last));
            fprintf('Iter %5d:  K = %4d (# change = %d)\n', ...
                t, numel(M), ch);             
        end
        
    end
    
    t = t + 1;
    
end


if displevel >= 1
    fprintf('Affinity propagation: # iters = %d,  converged = %d\n', ...
        t, converged);
end


% Extract results

M = get_centers(A, R).';
if ~isempty(M)
    K = numel(M);
    [~, L] = max(S(:,M), [], 2);
    L(M) = 1:K;
    L = L.';
else
    L = [];
end


%% core updating functions


function R = compute_r(S, A)

n = size(S, 1);
AS = A + S;
[Y, I] = max(AS, [], 2);

u = (1:n)' + (I - 1) * n;

AS(u) = -inf;
Y2 = max(AS, [], 2);

R = bsxfun(@minus, S, Y);
R(u) = S(u) - Y2;

function A = compute_a(R)

n = size(R, 1);
Rp = max(R, 0);

didx = 1 : (n+1) : n^2;
Rp(didx) = R(didx);

A = bsxfun(@minus, sum(Rp, 1), Rp);
dA = A(didx);
A = min(A, 0);
A(didx) = dA;

function M = get_centers(A, R)

M = find(diag(A) + diag(R) > 0);


function L = get_ass(S, M)

if ~isempty(M)
    K = numel(M);
    [~, c] = max(S(:,M), [], 2);
    c(M) = 1:K;
    L = M(c);
else
    n = size(S, 1);
    L = zeros(n, 1);
end



%% Option setting function

function [maxiter, lrate, exam_intv, displevel] = check_opts(nvlist)

maxiter = 500;
lrate = 0.5;
exam_intv = 10;
displevel = 0;

if ~isempty(nvlist)
    
    onames = nvlist(1:2:end);
    ovals = nvlist(2:2:end);
    
    for i = 1 : numel(onames)
        cn = onames{i};
        cv = ovals{i};
        
        switch lower(cn)
            case 'maxiter'
                if ~(isnumeric(cv) && isscalar(cv) && cv >= 1)
                    error('afpclus:invalidarg', ...
                        'MaxIter should be a positive integer.');
                end
                maxiter = cv;
                
            case 'learnrate'
                if ~(isfloat(cv) && isreal(cv) && isscalar(cv) && ...
                        cv > 0 && cv <= 1)
                    error('afpclus:invalidarg', ...
                        'LearnRate should be a real scalar in (0, 1].');
                end
                lrate = cv;
                
            case 'examintv'
                if ~(isnumeric(cv) && isscalar(cv) && cv == fix(cv) && cv >= 1)
                    error('afpclus:invalidarg', ...
                        'ExamIntv should be a positive integer.');
                end
                exam_intv = cv;
                
            case 'display'
                if ~ischar(cv)
                    error('afpclus:invalidarg', ...
                        'Display should be a char string.');
                end
                
                switch cv
                    case 'off'
                        displevel = 0;
                    case 'final'
                        displevel = 1;
                    case 'phase'
                        displevel = 2;
                    otherwise
                        error('afpclus:invalidarg', ...
                            'The value of Display is invalid.');
                end                
                
            otherwise
                error('afpclus:invalidarg', 'Unknown option name %s', cn);
        end
    
    end
    
end




