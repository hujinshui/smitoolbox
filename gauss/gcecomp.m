classdef gcecomp < gcebase
    % % The class for estimating generic complete covariance matrix
    %
    
    % Created by Dahua Lin, on April 14, 2010
    % Modified by Dahua Lin, on April 19, 2010
    %
    
    
    methods
                
        function C = estcov(obj, X, nw, mu) %#ok<MANU>
            % Implementation of the estimation algorithm
                                   
            d = size(X, 1);
            
            if isempty(nw) || size(nw, 1) == 1
                if isempty(nw)
                    EX2 = (X * X') * (1 / size(X,2));
                else
                    if d == 1
                        EX2 = (X .* nw) * X';
                    else
                        EX2 = bsxfun(@times, X, nw) * X';
                    end
                end
                
                if ~isequal(mu, 0)
                    S = EX2 - mu * mu';
                else
                    S = EX2;
                end
                S = 0.5 * (S + S');                
                
            else
                m = size(nw, 1);
                
                if d == 1
                    EX2 = (X.^2) * nw';
                    
                    if ~isequal(mu, 0)
                        S = EX2 - mu .^ 2;
                    else
                        S = EX2;
                    end
                    
                    S = reshape(S, [1 1 m]);
                else
                    S = zeros(d, d, m, class(X(1) * nw(1)));
                    zmu = isequal(mu, 0);
                    
                    for k = 1 : m
                        cS = bsxfun(@times, X, nw(k,:)) * X';
                        if ~zmu
                            cmu = mu(:,k);
                            cS = cS - cmu * cmu';
                            S(:,:,k) = 0.5 * (cS + cS');
                        end
                    end                    
                end                    
            end
            
            C = pdmat_gen(S);
        end
        
        
        function C = estcov_tied(obj, X, wi, wk, mu) %#ok<MANU>
            % Implementation of the estimation algorithm with tied covariance
                        
            Ex2 = bsxfun(@times, X, wi) * X';
            Eu2 = bsxfun(@times, mu, wk) * mu';
            
            S = Ex2 - Eu2;
            S = 0.5 * (S + S');
            
            C = pdmat_gen(S);
        end
        
    end
end