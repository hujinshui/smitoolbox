function v = logcrp(alpha, M, op)
% LOGCRP evaluates log-probability of a partition under Chinese Restaurant Process
%
%   Given a partition with K clusters, in which the k-th cluster has
%   m(k) elements, the log-probability of m is 
%
%       log p(m) = sum_{i=1}^K ( log(alpha) + log(gamma(m(k)-1)) )
%                - sum_{i=0}^{n-1} log(alpha + i)
%
%   In particular, we call the part of the formula in the first line
%   to be the "upper part", and the remaining part is called "lower part".
%   Note that the "lower part" is fixed when the total number of elements
%   is fixed, and thus sometimes we are only interested in the upper part.
%
%   v = LOGCRP(alpha, m);
%       evaluates the log-probability of a partition given by m under
%       a Chinese restaurant process of concentration parameter alpha.
%
%       m represents a partition as follows: there are numel(m) clusters,
%       the k-th cluster has m(k) elements. 
%
%   v = LOGCRP(alpha, m, 'u');
%       evaluates the upper part of the log-probability.
%
%   V = LOGCRP(alpha, M);
%   V = LOGCRP(alpha, M, 'u');
%
%       evaluates the log-probabilities of a collection of partitions
%       in M. Here M is a cell array, and each cell contains a partition.
%
%       The size of V is the same as the size of M.
%

% Created by Dahua Lin, on Sep 17, 2011
%

%% Verify input arguments

if ~(isfloat(alpha) && isscalar(alpha) && isreal(alpha) && alpha > 0)
    error('logcrp:invalidarg', 'alpha should be a positive real number.');
end

if isnumeric(M)
    m = M;
    if ~(isfloat(m) && isvector(m) && isreal(m))
        error('logcrp:invalidarg', 'Each m should be a real vector.');
    end
    n = 1;
elseif iscell(M)
    n = numel(M);
    for i = 1 : n
        m = M{i};
        if ~(isfloat(m) && isvector(m) && isreal(m))
            error('logcrp:invalidarg', 'Each m should be a real vector.');
        end
    end
else
    error('logcrp:invalidarg', ...
        'The 2nd argument should be either a real vector or a cell array.');
end
   
cl = true;
if nargin >= 3
    if ~(ischar(op) && strcmpi(op, 'u'))
        error('logcrp:invalidarg', 'The 3rd argument is invalid.');
    end
    cl = false;
end


%% Main

if n == 1
    v = eval_v(alpha, m, cl);
else
    v = zeros(size(M));
    for i = 1 : n
        v(i) = eval_v(alpha, M{i}, cl);
    end
end


%% Evaluation

function v = eval_v(alpha, m, cl)

K = numel(m);
v = K * log(alpha) + sum(gammaln(m));

if cl
    N = sum(m);
    vl = sum(log( alpha + (0:N-1) ));
    v = v - vl;
end







