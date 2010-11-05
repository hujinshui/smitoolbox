function g = parcorrgraph(J)
% construct partial correlation graph
%
%   g = parcorrgraph(J);
%       constructs a partial correlation graph given information matrix J.
%
%       In the ouput, g is a gr_adjlist struct. There is an edge between
%       i and j, when J(i,j) ~= 0, and the edge weight is 
%       -J(i,j) / sqrt(J(i,i) * J(j,j));
%

% Created by Dahua Lin, on Nov 5, 2010
%

%% verify input arguments

n = size(J, 1);
if ~(isfloat(J) && ndims(J) == 2 && n == size(J, 2))
    error('parcorrgraph:invalidarg', 'J should be a numeric matrix of size n x n.');
end

%% make graph

dv = full(diag(J));

if ~all(dv == 1)
    sca = 1 ./ sqrt(dv);
    J = bsxfun(@times, bsxfun(@times, J, sca), sca.');
end

[s, t, w] = find(J);
se = s < t;
s = s(se);
t = t(se);
w = - w(se);

g = gr_adjlist('u', n, s, t, w);



