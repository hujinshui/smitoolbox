classdef genmodel_base
    % The base of the classes that implement simple generative model
    %
    %   The model formulation:
    %       params ~ prior
    %       observations ~ params
    %
    %   This is an abstract class that defines the functions to be
    %   implemented by derived classes.
    %
    
    methods(Abstract)
        
        tf = is_supported_optype(model, optype);
        % Tests whether a particular operation type is supported
        %
        %   model.is_supported_optype('sample');
        %       returns whether (MCMC) sampling is supported.
        %
        %   model.is_supported_optype('atom');
        %       returns whether non-parameteric atom sampling is supported
        %
        %   model.is_supported_optype('varinfer');
        %       returns whether variational inference is supported.
        %
        
        n = get_num_observations(model, obs);
        % Gets the number of samples in the observation set
        %
        %   n = model.get_num_observations(obs);
        %       verifies the validity of obs as an observation set.
        %       If obs is valid, returns the number of samples in obs,
        %       otherwise, raise an error.
        %
        
        tf = is_valid_prior(model, pri);
        % Tests whether a given prior is a valid prior model
        %
        %   tf = model.is_valid_prior(pri);
        %       tests whether pri is a valid prior model of the
        %       parameters.
        %
        
        params = posterior_params(model, pri, obs, Z, optype);
        % Estimates/samples the parameters conditoned on given observations
        %
        %   atom = model.posterior_params(pri, obs, I, 'atom');
        %       draws an atom given a subset of data points selected by I.
        %
        %   params = model.posterior_params(pri, obs, Z, 'sample');
        %       draws K parameters given grouped/weighted data set.
        %
        %   params = model.posterior_params(pri, obs, Z, 'varinfer');
        %       estimates/infer K parameters given grouped/weighted data
        %       set.
        %
        %       Input Arguments:
        %       - pri:      the prior of the parameters
        %       - obs:      the observation data set
        %       - I:        an index vector.
        %       - Z:        It can be either a cell array of grouped
        %                   indices, or a K x n matrix of weights.
        %
        %       Output Arguments:
        %       - atom:     a single parameter.
        %
        %       - params:   The struct/object that captures the single or
        %                   multiple estimated parameters.
        %
        %       Note: atom and params can be in different forms.
        %
        
        Lliks = evaluate_logliks(model, params, obs);
        % Evaluates the logarithm of likelihood of samples
        %
        %   Lliks = model.evaluate_logliks(params, obs);
        %   Liiks = model.evaluate_logliks(atom, obs);
        %
        %       evaluates the log-likelihood values of the observations
        %       in obs with respect to the given set of parameters or
        %       the given atom.
        %
        %       If params has K parameters, Lliks is a K x n matrix.
        %
        %   Lliks = model.evaluate_logliks(pri, obs);
        %
        %       Under atom-mode, this evaluates the log marginal 
        %       likelihood values with respect to the given prior.    
        %
        %   Lliks = model.evaluate_logliks(.., obs, I);
        %
        %       Evaluates the log-likelihood values for a set of 
        %       observations selected by I.
        %
        
        Lpri = evaluate_logpri(model, pri, params);
        % Evaluate the log-prior of a given set of parameters
        %
        %   Lpri = model.evaluate_logpri(pri, params);
        %   Lpri = model.evaluate_logpri(pri, atom);
        %
        %       evaluates the log-prior of the given set of parameters
        %       or the given atom.
        %
        %       If params has K parameters, Lpri is a 1 x K vector.
        %
        
    end
    
end
