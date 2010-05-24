classdef gaussd_mle < smodelest
    % The Maximum likelihood estimator of Gaussian distribution
    %
    
    % Created by Dahua Lin, on Apr 14, 2010
    % Modified by Dahua Lin, on Apr 19, 2010
    %
               
    properties(GetAccess='public', SetAccess='private')
        cov_est;           % the covariance estimator
    end
      
    properties
        tie_cov = false;    % whether to tie covariance in estimation
        zero_mean = false;  % whether all mean vectors are zeros
        make_cp = false;    % whether to ouput object of class gaussd_cp
    end
                
    methods
        
        function obj = gaussd_mle(gce, varargin)
            % Constructs a Gaussian Maximum Likelihood estimator
            %
            %   obj = gaussd_mle(gce);
            %       constructs a Gaussian maximum likelihood estimator.
            %       Here, gce is the Gaussian covariance matrix estimator
            %       which is of class gmebase or its derived class.            
            %
            %   obj = gaussd_mle(gce, ...);
            %       One can further specify options to change the 
            %       default behavior of the estimator.
            %
            %       Each option is a string:
            %       - 'tie_cov':    tie the estimation of covariance
            %       - 'zero_mean':  assume zero mean vectors
            %       - 'make_cp':    output gaussd_cp object
            %
            
            if ~isa(gce, 'gcebase')
                error('gaussd_mle:invalidarg', ...
                    'gce should be an object of class gcebase.');                                
            end
            
            obj.cov_est = gce;     
            
            if ~isempty(varargin)                
                for i = 1 : length(varargin)
                    op = varargin{i};
                    
                    switch op
                        case 'tie_cov'
                            obj.tie_cov = true;
                        case 'zero_mean'
                            obj.zero_mean = true;
                        case 'make_cp'
                            obj.make_cp = true;
                        otherwise
                            error('gaussd_mle:invalidarg', ...
                                'Invalid option %s', op);
                    end                                        
                end                
            end
            
        end
        
        
        function Pb = accept(obj, X)
            % Checks the input observation and forms an estimation problem
            %                        
                        
            Pb = gaussd_mle_pb(obj, X);        
        end
                                       
    end    
            
end

