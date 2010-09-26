function d = clusvid(L1, L2)
% Compute the variation information distance between clusters
%
%   d = clusvid(L1, L2);
%       computes the variation information distance between two
%       clusters, which respectively label the same set of samples
%       with labels L1 and L2.
%
%       In the input L1 and L2 should be vectors of the same size.
%

%   History
%   -------
%       - Created by Dahua Lin, on May 26, 2010
%

%% verify input arguments

if ~(isvector(L1) && isvector(L2) && isnumeric(L1) && isnumeric(L2) && ...
        isequal(size(L1), size(L2)))
    error('clusvid:invalidarg', ...
        'L1 and L2 should be numeric vectors of the same size.');
end

%% main

b1 = min(L1) - 1;
b2 = min(L2) - 1;

if b1 ~= 0
    L1 = L1 - b1;
end

if b2 ~= 0
    L2 = L2 - b2;
end

m1 = max(L1);
m2 = max(L2);
n = numel(L1);

p1 = intcount([1, m1], L1) / n;
p2 = intcount([1, m2], L2) / n;
P12 = confusmat([m1, m2], L1, L2, 'nrm');

ev1 = calc_ent(p1);
ev2 = calc_ent(p2);
ev12 = calc_ent(P12);

d = 2 * ev12 - (ev1 + ev2);



%% subfunction

function ev = calc_ent(p)
% A sub-function to calculate entropy

p = p(p > 0);
ev = - sum(p .* log(p));

