function dists = raddiff(X1, X2)
%RADIANDIFF Computes the radian differences between corresponding vectors
%
%   dists = radiandiff(X1, X2);
%       computes the radian differences between corresponding vectors in X1
%       and X2.
%
%       The radian difference is the angle between two vectors in unit of
%       radian. In mathematics, it can be defined as arccos of normalized
%       correlation, which ranges from 0 (when correlation is 1) to pi
%       (when correlation is -1). If the radian difference between two
%       vectors is pi / 2, they are orthogonal to each other.
%
%       X1 and X2 should be matrices of the same size. Let their size be 
%       d x n, then dists will be a 1 x n vector, with dists(i) being the
%       radian difference betwen X1(:, i) and X2(:, i).
%

%   History
%       - Created by Dahua Lin, on Jun 3, 2008
%       - Modified by Dahua Lin, on Jul 22, 2010
%           - based on nrmdot
%       - Modified by Dahua Lin, on Aug 2, 2010
%           - rename: radiandiff -> raddiff
%

%% main

nds = nrmdot(X1, X2);

nds(nds > 1) = 1;
nds(nds < -1) = -1;

dists = acos(nds);


