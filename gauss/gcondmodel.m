classdef gcondmodel < modelest
    % The abstract base class for all Gaussian conditional models
    %
    %   A Gaussian conditional model is a likelihood model that
    %   is conditioned on Gaussian-distributed variables.
    %
    %   
    
    % Created by Dahua Lin, on April 14, 2010
    % Updated by Dahua Lin, on April 19, 2010
    %   - the estimate functionalities are discarded, which are moved
    %     to the new class gcondmle.
    %
    
    properties(Constant)
        est_type = 'mle';
    end    
        
    properties(GetAccess='public', SetAccess='private')
        cdim;       % the dimension of the variables being conditioned
    end    
    
    methods
        function obj = gcondmodel(cdim)
            % Constructs a Gaussian conditional model
            %
            %   obj = gcondmodel(cdim);
            %       constructs a Gaussian conditional model in which
            %       the dimension of the variables being conditioned
            %       is cdim.            
            %
            
            obj.cdim = cdim;
        end
        
        function R = estimate(obj, X, W)
            % Perform MLE estimation of the hidden variable
            %
            %   R = obj.estimate(X, W);
            %       performs maximum a likelihood estimation of the
            %       hidden variables that the observations are conditioned
            %       upon.
            %
            %       It invokes the the conj_update to compute U1 and U2
            %       and then compute U2 \ U1.
            %
            
            if nargin < 3
                W = [];
            end            
            [U1, U2] = conj_update(obj, X, W);
            R = cldiv(U2, U1);
        end
                                
    end    
    
    methods(Abstract)
        
        [U1, U2] = conj_update(obj, X, W);
        % Compute the conjugate update of the canonical parameters of 
        % the Gaussian prior based on input observations
        %
        % In the input, X gives the observations, whose form are 
        % model-dependent. W gives the weights of each observation.
        %
        % In the output, U1 and U2 are respectively the terms for 
        % updating theta1 and theta2.
        %
        
        L = eval_loglik(obj, M, X);
        % Evaluate the log likelihood
        %
        %   L = obj.eval_loglik(M, X);
        %       evaluates the log likelihood of the samples in X
        %       with respect to the models with the hidden variables
        %       given by columns of M.
        %
        %       If the size of M is d x m, and that of X is d x n,
        %       then L is of size m x n.
        %
                
    end    
    
end