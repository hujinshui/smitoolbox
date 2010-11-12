function X = ddsample(P, n, rstream)
% Draw samples from discrete distribution(s)
%
%   X = ddsample(P, n);
%   X = ddsample(P, n, rstream);
%       draw n samples from each of the discrete distributions in P.
%
%       Let P be an K x m matrix, then P represents m distinct discrete
%       distributions over K classes. 
%
%       The output X will be an n x m matrix, where X(:,j) contains the
%       samples drawn from the distribution P(:,j).
%
%       One can also specify the random number stream used for generating
%       random numbers.
%

%   Created by Dahua Lin, on Nov 7, 2010
%

%% verify input arguments

if ~(isfloat(P) && ~issparse(P) && ndims(P) == 2)
    error('ddsample:invalidarg', ...
        'P should be a non-sparse numeric matrix.');
end

if ~(isnumeric(n) && isscalar(n) && n >= 1)
    error('ddsample:invalidarg', ...
        'n should be an integer scalar.');
end
n = double(n);

if nargin < 3
    rstream = RandStream.getDefaultStream;
end

%% main

F = cumsum(P, 1);

V = rand(rstream, n, size(P,2));
X = ddsample_cimp(F, V);

