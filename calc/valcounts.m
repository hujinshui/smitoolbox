function [C, U] = valcounts(A, U, op)
%VALCOUNTS Counts the occurrences of distinct values
%
%   [C, U] = valcounts(A);
%       counts the number of occurrences of distinct values in A. 
%       U is an array of unique values, C is an numeric array of the
%       same size as U, which gives the corresponding counts.
%
%   [C, U] = valcounts(A, [], 'rows');
%       counts the number of occurrences of distinct rows in A.
%
%   C = valcounts(A, U);
%       counts the number of the occurrences of the values given in U.
%
%   C = valcounts(A, U, 'rows');
%       counts the number of the occurrences of rows given in U.
%

%   History
%       - Created by Dahua Lin, on May 31, 2008
%

%% parse and verify input arguments

if nargin < 2
    U = [];
end

if nargin >= 3
    assert(isequal(op, 'rows'), 'valcounts:invalidarg', ...
        'Invalid input arguments to valcounts.');
    
    for_rows = true;
else
    for_rows = false;
end

%% main

if isempty(U)
   
    if ~for_rows 
        [U, ~, J] = unique(A);
    else
        [U, ~, J] = unique(A, 'rows');
    end
    
    Js = sort(J);
    [sp, ep] = valueseg(Js);    
    C = ep - sp + 1;
    
else
    
    if ~for_rows
        [~, J] = ismember(A, U);
        C = zeros(size(U));
    else
        [~, J] = ismember(A, U, 'rows');
        C = zeros(size(U, 1), 1);
    end
    
    Js = sort(J(J > 0));
    [sp, ep, jv] = valueseg(Js);
    C(jv) = ep - sp + 1;
    
end


