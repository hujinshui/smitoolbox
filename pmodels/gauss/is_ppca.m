function tf = is_ppca(M)
% Tests whether the input argument is a PPCA model struct
%
%   tf = is_ppca(M);
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% main

tf = isstruct(M) && numel(M) == 1 && isfield(M, 'tag') && ...
    strcmp(M.tag, 'ppca');

