classdef multd_est < smodelest
    % The class for maximum a posteriori estimation of multi-class 
    % discrete distribution
    %
    % The estimation can use Dirichlet distribution as prior. The 
    % observation is a q-map, i.e. a soft assignment of each sample 
    % to different classes
    %
    
    % Created by Dahua Lin, on April 22, 2010
    %
    
        
    methods
        
        function obj = multd_est(pri)
            % Constructs a binary distribution estimator
            %
            %   obj = binaryd_map();
            %       constructs a maximum likelihood estimator of binary 
            %       distribution.            
            %
            %   obj = binaryd_map(pri);
            %       constructs a maximum a posteriori estimator of 
            %       binary distribution.
            %
            
            if nargin >= 1 && ~isempty(pri)
                if ~(isa(pri, 'dirichletd') && pri.nmodels == 1)
                    error('multd_est:invalidarg', ...
                        'pri should be an object of class dirichled.');
                end
                
                obj.prior = pri;
            end                                            
        end
        
        
        function Pb = accept(obj, Q)
            % Checks the input observation and forms an estimation problem
            %                        
                        
            Pb = multd_est_pb(obj, Q);
        end

    end    
    
end