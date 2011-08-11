function sm = sm_select(sm, i)
% select a subset of matrices from an symat object
%
%   sm_r = sm_select(sm, i);
%       selects a subset (specified via i, an index or a vector of index)
%       of matrices from sm to form a new symat object, returned as sm_r.
%

% Created by Dahua Lin, on Aug 11, 2011
%

switch sm.ty
    case {'s', 'd'}
        sm.v = sm.v(:, i);
        sm.n = size(sm.v, 2);
    case 'f'
        sm.v = sm.v(:,:,i);
        sm.n = size(sm.v, 3);
end

