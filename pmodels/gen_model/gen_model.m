classdef gen_model
    % The base class for generative model
    %
    % In general, the math formulation of the model can be written as
    %
    %   params ~ prior distribution ( hyper-params )
    %   output-vars ~ distribution ( params, input-vars, labels )
    %
    %   Specifically, each output variable is generated conditioned on
    %   the corresponding input variables (or none), based on the
    %   model whose parameters are selected by corresponding label.
    %
    
    % Created by Dahua Lin, on Aug 26, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')                
        input_vars_info;        % the information of input variables
        output_vars_info;       % the information of output variables
        params_info;            % the information of parameters
        hyper_params_info;      % the information of hyper-parameters        
    end
    
    properties(Dependent)        
        num_input_vars;     % the number of input variables
        num_output_vars;    % the number of output variables
        num_params;         % the number of parameters
        num_hyper_params;   % the number of hyper parameters        
    end
        
    methods
        function n = get.num_input_vars(model)
            n = numel(model.input_vars_info);
        end
        
        function n = get.num_output_vars(model)
            n = numel(model.output_vars_info);
        end
        
        function n = get.num_params(model)
            n = numel(model.params_info);
        end
        
        function n = get.num_hyper_params(model)
            n = numel(model.hyper_params_info);
        end
    end    
        
    
    %% Constructors
    
    methods
        
        function model = gen_model(ivinfo, ovinfo, painfo, hpinfo)
            % Constructs the base of a generative model
            %
            %   model = gen_model(prinfo, painfo, hpinfo);
            %
            %   Input arguments:
            %   - ivinfo:   the input variable information struct
            %   - ovinfo:   the output variable information struct
            %   - painfo:   the param information struct
            %   - hpinfo:   the hyper parameter information struct
            %
            
            model.input_vars_info = ivinfo;
            model.output_vars_info = ovinfo;
            model.params_info = painfo;
            model.hyper_params_info = hpinfo;
        end
        
    end
    
       
    %% evaluation methods
    
    methods(Abstract)
    
        lpri = logpri(model, hparam, params, pid, hmap)
        % Evaluates the log-prior of given parameters
        %
        %   lpri = model.logpri(hparams, params, pid);
        %
        %       evaluates the log-prior of the given parameters (params),
        %       whose id is given by pid, with respect to the prior model 
        %       whose parameters are given by the hyper parameters 
        %       (hparams).,        
        %
        %   lpri = model.logpri(hparams, params, pid, hmap);
        %       
        %       evaluates the log-prior of the given parameters (params)
        %       with respect to multiple prior models whose parameters 
        %       are from hparams. 
        %
        %       Specifically, the log-prior of the k-th parameter should
        %       be evaluated based on the hmap(k)-th hyper-parameter.
        %
        %   For both syntax, the output lpri is a 1 x m vector, where
        %   m is the number of parameters given in params.
        %
        
        Llik = loglik(model, params, invars, outvars)
        % Evaluates the generative log-likelihood 
        %
        %   Llik = model.loglik(params, invars, outvars);
        %
        %       evaluates the log-likelihood of the observed variables
        %       given by invars (input) & outvars (output), with respect
        %       to the models whose parameters are given by params.
        %
        %       The evaluation is done in a pairwise manner. Specifically, 
        %       suppose there are n observations, and m sets of model
        %       parameters, then Llik is a matrix of size m x n, such 
        %       that Llik(k, i) is the likelihood of the i-th observation
        %       with respect to the k-th parameter.
        %               
        
    end
        
    
    %% sampling & inference methods
    
    methods(Abstract)
               
        outvars = generate(model, invars, params, n, g);
        % Generates output variables given inputs, parameters, and labels.
        %
        %   outvars = model.generate(invars, params, n);
        %
        %       generates n samples of the output variables from a single
        %       model whose parameter is given by params.
        %
        %   outvars = model.generate(invars, params, n, g);
        %
        %       generates n samples from multi-models whose parameters 
        %       are given by params, based on the grouping given in g, 
        %       a cell array of index vectors. 
        %
        %       In particular, the samples with indices in g{k} should
        %       be generated by the k-th model. Here, n is the total
        %       number of samples.
        %
        %       outvars is a cell array of size 1 x num_output_vars.
        %       outvars{id} is the set of outputs whose identifier is id.
        %
        
        theta = solve_params(model, X, Z, alpha, hmap);
        % Solves optimal parameters
        %
        %   params = model.solve_params(X);        
        %   params = model.solve_params(X, [], alpha);
        %
        %       Performs MAP estimation of the parameter based on 
        %       a given set of samples, each with weight 1. 
        %
        %       Alpha is a cell array comprised of all hyper parameters.
        %       When num_hyper_params == 0, Alpha should be empty or 
        %       omitted, and this function should perform MLE estimation.
        %
        %   Theta = model.mapest_params(X, W, alpha);
        %
        %       Performs MAP estimation based on weighted samples.
        %
        %       W should be a K x n matrix, where W(k, i) is the
        %       contribution weight of the i-th sample to the k-th model.
        %       In output, Theta is comprised of K parameters.
        %
        %   Theta = model.mapest_params(X, g, alpha);
        %
        %       Performs MAP estimation based on grouped samples.
        %       Here, g is a cell array of index vectors, and the k-th
        %       model parameter should be estimated based on the samples
        %       whose index is in g{k}.
        %
        %   Theta = model.mapest_params(X, .., Alpha, hmap);
        %
        %       Performs MAP estimation with hyper-map. Here, hmap is
        %       an index vector of length K, and Alpha is a cell matrix,
        %       of which each column gives a set of hyper-params.
        %
        %       The estimation of the k-th model parameter should be
        %       based on hmap(k)-th hyper-parameter(s) contained in
        %       Alpha. 
        %
        
        Theta = sample_params(model, X, Z, Alpha, hmap);
        % Samples from the posterior distribution of parameters
        %
        %   theta = model.sample_params(X, [], alpha);
        %
        %       samples a parameter theta, given a set of data X and
        %       the hyper parameter alpha.
        %
        %   theta = model.sample_params(X, [], alpha, n);
        %   theta = model.sample_params(X, w, alpha, n);
        %       
        %       samples n parameters from a given set of data X,
        %       which may or may not be weighted, and a single
        %       hyper-parameter alpha. 
        %
        %   Theta = model.sample_params(X, W, alpha);
        %   Theta = model.sample_params(X, g, alpha);
        %   Theta = model.sample_params(X, .., Alpha, hmap);
        %
        %       samples K parameters (packed in Theta as output), 
        %       given a weighted/grouped set of data X.
        %       hyper-map can also be used here.                
        %
        
    end
    
       
end



