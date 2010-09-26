function [Y, evs] = lle_embed(A, d)
% Solves the embedded coordinates with an LLE graph
%
%   [Y, evs] = lle_embed(A, d);
%       solves the embedded coordinates with respect to the given LLE 
%       coefficient graph A. Here, A(i,:) is the coeffcient vector for
%       constructing the i-th sample.
%
%       Here, d is the dimension of the embedded space.
%
%       In the output, Y is the coordinate matrix of size d x n, and
%       evs is a d x 1 vector comprised of the corresponding eigenvalues.

% History
% -------
%   - Created by Dahua Lin, on Nov 26, 2009
%

%% parse and verify input arguments

assert(isfloat(A) && isreal(A) && ndims(A) == 2 && size(A,1) == size(A,2), ...
    'lle_embed:invalidarg', 'A should be a real value square matrix.');

n = size(A, 1);

assert(isnumeric(d) && isscalar(d) && d == fix(d) && d >= 1 && d < n, ...
    'lle_embed:invalidarg', 'd should be a integer scalar in [1, n).');

%% main

ImA = speye(n) - A;
M = ImA' * ImA;

opts.issym = 1;
opts.disp = 0;
[V, D] = eigs(M, d+1, 'SM', opts);

Y = V';
evs = diag(D);

[evs, si] = sort(evs);
evs = evs(2:end);
Y = Y(si(2:end), :);

