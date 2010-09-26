function V = measure_mtimes_imp(P, X)
%MEASURE_MTIMES_IMP implements P' * X with the rule 0 * inf = 0
%
%   P and X should be float matrices with size(P,1) == size(X,1)
%
% Created by Dahua Lin, on Jun 5, 2008
%

%% main

assert(size(P,1) == size(X,1));

all_positive_p = all(P > 0, 1);

if all(all_positive_p > 0)
    
    V = P' * X;
    
else    
    m = size(P, 2);
    n = size(X, 2);
       
    V = zeros(m, n, class(P));    
    
    if any(all_positive_p)
        
        % treat those p with all positive values in batch
        
        V(all_positive_p, :) = P(:, all_positive_p)' * X;
        
        % treat those p with some zero values specially 
        
        V(~all_positive_p, :) = measure_mtimes_cimp(P(:, ~all_positive_p), X)';
                
    else
        
        V = measure_mtimes_cimp(P, X)';
    
    end
    
end


