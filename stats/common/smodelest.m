classdef smodelest
    % The abstract base class for statistical model estimator
    %
    
    % Created by Dahua Lin, on April 19
    %
    
    properties(GetAccess='public', SetAccess='protected')
        prior;
    end
            
    methods
        
        function M = estimate(obj, X, W)
            % Estimates the model(s) from observations
            %
            %   M = obj.estimate(X);
            %   M = obj.estimate(X, W);
            %       The function estimates the model(s) from (weighted)
            %       observations
            %
            %       The observed samples are given by X, while the
            %       weights of the samples are given by W.
            %       If there are n samples, then W can be in either of
            %       the following forms:
            %       - W can be empty or omitte. In this case, each sample 
            %         is with a unit weight.
            %       - W can be a row vector that gives the weights of
            %         all samples respectively.
            %       - W can be a m x n matrix where W(k, :) gives the
            %         weights for estimating the k-th model. In this case,
            %         m models are estimated, and M should be an object
            %         containing m models respectively estimated based on
            %         different sets of weights.
            %
            %       This function is just a short-cut of
            %
            %           obj.accept(X).estimate(W);
            %
            
            if nargin < 3
                W = [];
            end            
            M = obj.accept(X).estimate(W);
        end
    end
    
    
    methods
                        
        Pb = accept(obj, X);
        % Checks the input observation and forms an estimation problem
        %
        %   Pb = accept(obj, X);
        %       The function checks the integrity of the input 
        %       observations (X), and forms an estimation problem,
        %       which is an object of class smodelestpb.
        %
        %       Pre-computation that relies only on X but not on W
        %       can be performed at this function. In this way,
        %       when only the weight changes, the time for computing
        %       these quantities can be saved.
        %
        
    end
end