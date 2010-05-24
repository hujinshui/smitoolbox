classdef gcebase 
    % The base class for Gaussian covariance estimation
    %
    
    % Created by Dahua Lin, on Apr 14, 2010
    % Modified by Dahua Lin, on Apr 19, 2010
    %
        
    properties(GetAccess='public', SetAccess='protected')
        dim;    % if specified, it means that it can only estimate the 
                % the models of specific dimension. 
    end
           
    methods(Abstract)
        
        C = estcov(obj, X, nw, mu);
        % Estimates covariance using zero-mean observations
        %   
        %   C = obj.estcov(X, nw, mu);
        %
        % Input:
        %   - X:    The sample matrix, size d x n. 
        %   - nw:   normalized weights 
        %           empty (means each sample is with weight 1/n), or
        %           m x n, (each row sums to 1)
        %   - mu:   the pre-computed mean.
        %
        % Output:
        %   - covariance matrices as a pdmat object
        %
        
        C = estcov_tied(obj, X, wi, wk, mu);
        % Implementation of the estimation algorithm with tied covariance
        %
        %   C = obj.estcov_tied(X, wi, wk, mu);
        %
        % Input:
        %   - X:    sample matrix, size d x n. 
        %   - wi:   normalized aggregate per-sample weights [1 x n]
        %   - wk:   normalized aggregate per-model weights [1 x m]
        %   - mu:   the pre-computed mean.
        %
        % Output:
        %   - covariance matrix as a pdmat object
        %
        
    end
    
end

