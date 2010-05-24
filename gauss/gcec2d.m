classdef gcec2d < gcebase
    % % The class for estimating complete 2D covariance matrix
    %
    
    % Created by Dahua Lin, on April 14, 2010
    % Modified by Dahua Lin, on April 19, 2010
    %
    
    methods
        
        function obj = gcec2d()
            obj.dim = 2;
        end
                
        function C = estcov(obj, X, nw, mu) %#ok<MANU>
            % Implementation of the estimation algorithm
               
            x1 = X(1, :);
            x2 = X(2, :);
            
            x11 = x1 .^ 2;
            x12 = x1 .* x2;
            x22 = x2 .^ 2;
            
            if isempty(nw)                
                EX2 = [sum(x11); sum(x12); sum(x22)] / size(X, 2);
            else
                EX2 = [x11; x12; x22] * nw';
            end
               
            if ~isequal(mu, 0)
                u1 = mu(1, :);
                u2 = mu(2, :);
                
                v = EX2 - [u1 .^ 2; u1 .* u2; u2 .^ 2];
            else
                v = EX2;
            end
            
            C = pdmat_c2d(v);
        end
        
        
        function C = estcov_tied(obj, X, wi, wk, mu) %#ok<MANU>
            % Implementation of the estimation algorithm with tied covariance
               
            x1 = X(1, :);
            x2 = X(2, :);
            
            u1 = mu(1, :);
            u2 = mu(2, :);
                        
            v = [x1.^2; x1.*x2; x2.^2] * wi' - [u1.^2; u1.*u2; u2.^2] * wk';
            
            C = pdmat_c2d(v);
        end
        
    end
end