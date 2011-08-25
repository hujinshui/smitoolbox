function tf = is_pdmat(S)
% Tests whether an argument is a pdmat struct
%
%   tf = is_pdmat(S);
%       performs a quick check, and returns whether S is a pdmat struct.
%

% Created by Dahua Lin, on Aug 25, 2010
%

tf = isstruct(S) && isscalar(S) && isfield(S, 'tag') && ...
    ischar(S.tag) && strcmp(S.tag, 'pdmat');

