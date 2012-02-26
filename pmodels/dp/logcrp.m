function v = logcrp(alpha, M)
% LOGCRP evaluates log-probability of a partition under Chinese Restaurant Process
%
%   Given a partition with K clusters, in which the k-th cluster has
%   m(k) elements, the log-probability of m is 
%
%       p(m) = alpha^K * prod_{k=1}^K (m(k) - 1)! * 
%              Gamma(alpha) / Gamma(alpha + N).
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
%   V = LOGCRP(alpha, M);
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
    if ~(isfloat(M) && isvector(M) && isreal(M))
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
   

%% Main

if n == 1
    v = eval_v(alpha, M);
else
    v = zeros(size(M));
    for i = 1 : n
        v(i) = eval_v(alpha, M{i});
    end
end


%% Evaluation

function v = eval_v(alpha, m)

N = sum(m);
c = gammaln(alpha) - gammaln(alpha + N);

m = m(m > 0);
v = numel(m) * log(alpha) + sum(gammaln(m)) + c;




