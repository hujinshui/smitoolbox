classdef gcediag < gcebase
    % The class for estimating diagonal covariance matrix
    %
    
    % Created by Dahua Lin, on April 14, 2010
    % Modified by Dahua Lin, on April 19, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        iso;        % whether to restrict the covariance be isotropic
    end    
    
    methods
        
        function obj = gcediag(op)
            % constructs a Gaussian MLE estimator
            %
            %   obj = gmediag();
            %       constructs a Gaussian MLE estimator for estimating
            %       Gaussian distribution with diagonal covariance.
            %
            %   obj = gmediag('iso');
            %       constructs a Gaussian MLE estimator for estimating
            %       Gaussian distribution with isotropic covariance.
            %       
            
            if nargin >= 1
                if ~(ischar(op) && strcmp(op, 'iso'))
                    error('gcediag:invalidarg', ...
                        'The first argument can only be ''iso''.');
                end
                use_iso = true;
            else
                use_iso = false;
            end
               
            obj.iso = use_iso;
        end
        
        
        function C = estcov(obj, X, nw, mu)
            % Estimating covariance (respectively)
                                    
            if isempty(nw)
                EX2 = sum(X.^2, 2) / size(X, 2);
            else
                EX2 = (X.^2) * nw';
            end
               
            if ~isequal(mu, 0)
                dv = EX2 - mu.^2;
            else
                dv = EX2;
            end
            
            if obj.iso
                d = size(dv, 1);        
                if d == 1
                    C = pdmat_udiag(d, dv);
                else
                    C = pdmat_udiag(d, sum(dv, 1) / d);
                end
            else
                C = pdmat_diag(dv);
            end                            
        end
        
        
        function C = estcov_tied(obj, X, wi, wk, mu)
            % Estimating the tied covariance.
                        
            dv = (X.^2) * wi' - (mu.^2) * wk';
            
            if obj.iso
                d = size(dv, 1);
                if d == 1
                    C = pdmat_udiag(d, dv);
                else
                    C = pdmat_udiag(d, sum(dv, 1) / d);
                end
            else
                C = pdmat_diag(dv);
            end
        end
        
    end
end