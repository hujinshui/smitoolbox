function sm = sm_scale(sm, c)
% Scale the symmetric matrices in an symat struct
%
%   sm = sm_scale(sm, c);
%       multiples the symmetric matrices in sm with coefficient c.
%

% Created by Dahua Lin, on Aug 11, 2011
%


if ~(isfloat(c) && isscalar(c))
    error('sm_scale:invalidarg', 'The coefficient should be a numeric scalar.');
end

sm.v = sm.v * c;
