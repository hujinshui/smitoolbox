function [x, acc] = crpsim(alpha, n, pricount)
%CRPSIM Simulate a Chinese Restaurant Process (CRP)
%
%   x = CRPSIM(alpha, n);
%
%       simulates a Chinese Restaurant process to generate a sequence
%       of length n.
%
%   x = CRPSIM(alpha, n, pricount);
%
%       performs the simulation with prior count on the first K values.
%       K is numel(pricount).
%
%   [x, acc] = CRPSIM(alpha, n);
%   [x, acc] = CRPSIM(alpha, n, pricount);
%
%       Additionally returns the accumulated counts of each atom. 
%       If pricount is given, acc includes the prior counts.
%

% Created by Dahua Lin, on Sep 17, 2011
%

%% verify input arguments

if ~(isfloat(alpha) && isscalar(alpha) && isreal(alpha) && alpha > 0)
    error('crpsim:invalidarg', 'alpha should be a positive real number.');
end

if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
    error('crpsim:invalidarg', 'n should be a non-negative integer.');
end

if nargin < 3
    pricount = [];
else
    if ~(isnumeric(pricount) && isvector(pricount) && isreal(pricount))
        error('crpsim:invalidarg', 'pricount should be a real vector.');
    end
end

%% main

alpha = double(alpha);
if ~isempty(pricount) && ~isa(pricount, 'double')
    pricount = double(pricount); 
end

r = rand(1, n);
if nargout < 2
    x = crpsim_cimp(alpha, r, pricount);
else
    [x, acc] = crpsim_cimp(alpha, r, pricount);
end

