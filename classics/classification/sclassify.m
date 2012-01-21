function [L, v] = sclassify(D, op)
%SCLASSIFY Standard Nearest Neighbor Classification
%
%   L = SCLASSIFY(D, 'dis');
%   L = SCLASSIFY(D, 'sim');
%
%       Classifies each sample to a particular class that is closest to it.
%
%       Suppose there are K class and n testing samples. Then D should be
%       matrix of size K x n, and D(k, i) is the distance (2nd arg being
%       'dis') or similarity (2nd arg being 'sim') between the i-th sample
%       and the k-th class.
%
%       L is a vector of size 1 x n, whose values are integers in 
%       {1, ..., K}.
%
%
%   [L, v] = SCLASSIFY( ... );
%     
%       Additionally returns the distance/similarity value between each
%       sample to the selected class.
%
%

% Created by Dahua Lin, on Jun 6, 2010
% Modifued by Dahua Lin, on Jan 21, 2012
%

%% verify input

if ~(isnumeric(D) && ndims(D) == 2)
    error('sclassify:invalidarg', 'D should be a numeric matrix.');
end

if strcmp(op, 'dis')
    is_dis = true;
elseif strcmp(op, 'sim')
    is_dis = false;
else
    error('sclassify:invalidarg', ...
        'The second arg should be either ''dis'' or ''sim''.');
end
   
%% main

if is_dis    
    [v, L] = min(D, [], 1);
else
    [v, L] = max(D, [], 1);
end

    
    
