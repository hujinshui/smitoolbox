function tf = is_svm_dual_sol(S)
%IS_SVM_PROBLEM Whether the input argument is an SVM dual-solution struct
%
%   tf = is_svm_dual_sol(S);
%       returns whether S is a valid SVM problem struct.
%

% Created by Dahua Lin, on Jan 17, 2012
%

%% main

tf = isstruct(S) && numel(S) == 1 && isfield(S, 'tag') && ...
    strcmp(S.tag, 'svm-dual-sol');
