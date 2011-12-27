function [zty, K] = verify_Zarg(Z, n)
% Verifies the Z argument for sample weighting/grouping
%
%   [zty, K] = verify_Zarg(Z, n);
%
%       Here n is the total number of samples
%   
%       If Z is empty, then returns zty = 0 and K = 1.
%
%       If Z is a weight matrix, then returns zty = 1 and K = size(Z,1).
%
%       If Z is a cell array of index vectors, then returns zty = 2 and
%       K = numel(Z).
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% main

zty = -1;

if isempty(Z)
    zty = 0;
    K = 1;
    
elseif isnumeric(Z)
    if isfloat(Z) && isreal(Z) && ndims(Z) == 2 && size(Z,2) == n
        zty = 1;
        K = size(Z, 1);
    end
    
elseif iscell(Z)
    if isvector(Z) && all(cellfun(@isnumeric, Z))
        zty = 2;
        K = numel(Z);
    end
end

if zty < 0
    error('verify_Zarg:invalidarg', 'The Z argument is invalid.');
end


    
    