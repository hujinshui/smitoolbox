classdef binaryd_est < smodelest
    % The class for maximum a posteriori estimation of binary distribution
    %
    % The estimation can use beta distribution as prior. The observation
    % is a row vector that gives the soft assignment of each sample to
    % label 1.
    %
    
    % Created by Dahua Lin, on April 22, 2010
    %
    
        
    methods
        
        function obj = binaryd_est(pri)
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
                if ~(isa(pri, 'betad') && pri.nmodels == 1)
                    error('binaryd_est:invalidarg', ...
                        'pri should be an object of class betad.');
                end
                
                obj.prior = pri;
            end                                            
        end
        
        
        function Pb = accept(obj, Q)
            % Checks the input observation and forms an estimation problem
            %                        
                        
            Pb = binaryd_est_pb(obj, Q);
        end

    end    
    
end