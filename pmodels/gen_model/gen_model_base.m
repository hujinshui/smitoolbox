classdef gen_model_base
    % The base class of all generative models
    %
    % Generally, the model is formulated as follows
    %
    %   - a parameter is generated from a prior distribution 
    %     characterized by relevant hyper parameters
    %
    %   - a product sample is generated from the corresponding
    %     parameteric model. The corresponding parameters can
    %     be selected via a label. 
    %
    
    % Created by Dahua Lin, on Aug 28, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')        
        products_info;          % the information of products
        params_info;            % the information of parameters
        hyper_params_info;      % the information of hyper-parameters        
    end
    
    properties(Dependent)        
        num_products;       % the number of products
        num_params;         % the number of parameters
        num_hyper_params;   % the number of hyper parameters        
    end
        
    methods
        function n = get.num_products(model)
            n = numel(model.products_info);
        end
        
        function n = get.num_params(model)
            n = numel(model.params_info);
        end
        
        function n = get.num_hyper_params(model)
            n = numel(model.hyper_params_info);
        end
    end
    
    
    %% Constructor
    
    methods
        
        function model = gen_model_base(Pr, Pa, HPa)
            % Constructs the base of a generative model
            %
            %   model = gen_model_base(Pr, Pa, HPa);
            %       
            %   Input arguments:
            %       - Pr:   the info of products
            %       - Pa:   the info of parameters
            %       - HPa:  the info of hyper-parameters
            %
            %   Each of this argument is a struct array, which
            %   contains the following fields: name, type, and size.
            %
            
            model.products_info = Pr;
            model.params_info = Pa;
            model.hyper_params_info = HPa;
        end
        
    end
    
    %% Interfaces
    
    methods
        
        LLik = loglik(model, params, obs);
        % Evaluates the logarithm of likelihood values
        %
        %   LLik = loglik(model, params, obs);
        %
        %       Evaluates the logarithm of the likelihood values of
        %       observations with respect to models, in a pairwise
        %       manner.
        %
        %       Specifically, if there are m parameters in params, 
        %       and n observations in obs. The output LLik should be
        %       a matrix of size m x n, in which LLik(k, i) is the
        %       log-likelihood value of the i-th observation with
        %       respect to the k-th model parameter.
        %        
        
    end
        
    
end



