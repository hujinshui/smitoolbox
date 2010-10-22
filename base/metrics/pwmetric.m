function dists = pwmetric(metricfun, X1, X2)
%PWMETRIC Computes the pairwise metric matrix based on a supplied metric
%function
%
%   dists = pwmetric(metricfun, X1, X2);
%       computes the pairwise metric matrix based on metricfun.
%
%       The metricfun should be a function supporting the following syntax:
%
%           v = metricfun(Y1, Y2);
%
%       with Y1 and Y2 both being d x n matrices, v should be a 1 x n
%       vector, with v(i) = d(Y1(:, i), Y2(:, i));
%
%       This function uses metricfun to compute pairwise metrics, and
%       outputs dists, which has
%
%           dists(i, j) = d(X1(:, i), X2(:, i));
%
%       Suppose X1 and X2 respectively have n1 and n2 columns, then 
%       dists is an n1 x n2 matrix.
%
%   dists = pwmetric(metricfun, X);
%       computes pairwise metrics between the columns in X, which is
%       equivalent to pwmetric(metricfun, X, X);
%

%   History
%       - Created by Dahua Lin, on Jun 2, 2008
%       - Modified by Dahua Lin, on Jun 3, 2008
%           - automatic decide the class of output based on first batch
%             of results.
%

%% parse and verify input arguments

if nargin < 3
    X2 = X1;
end

assert(ndims(X1) == 2 && ndims(X2) == 2 && size(X1, 1) == size(X2, 1), ...
    'pwmetric:invalidarg', ...
    'X1 and X2 should be matrices with the same length along 1st dimension.');


%% main

n1 = size(X1, 2);
n2 = size(X2, 2);

if n1 < n2  
    
    y = X1(:, 1);
    r1 = metricfun(y(:, ones(1, n2)), X2);
    dists = zeros(n1, n2, class(r1));
    dists(1, :) = r1;
    
    for i = 2 : n1        
        y = X1(:, i);
        dists(i, :) = metricfun(y(:, ones(1, n2)), X2);        
    end    
    
else  
    
    y = X2(:, 1);
    r1 = metricfun(X1, y(:, ones(1, n1)))';
    dists = zeros(n1, n2, class(r1));
    dists(:, 1) = r1;
    
    for i = 2 : n2
        y = X2(:, i);
        dists(:, i) = metricfun(X1, y(:, ones(1, n1)))';
    end
    
end



