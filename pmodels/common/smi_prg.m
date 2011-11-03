classdef smi_prg
    % The base class for statistical modeling and inference program
    %
    %   To implement a probabilistic inference and estimation framework,
    %   one can write a class that derives from this, and use a driver
    %   function such as smi_mcmc or smi_varinfer to drive the inference
    %   and estimation procedure. In this way, the framework writer does
    %   not need to write his own procedural control code. 
    %   
   
    % Created by Dahua Lin, on Sep 3, 2011
    %
    
    methods(Abstract)
        
        [Sd, Sc] = initialize(prg, obs, S0, optype);
        % Initialize the states given observations and other fixed params
        %
        %   [Sd, Sc] = prg.initialize(obs, S0, optype);
        %       
        %       This method accepts an observation struct which contains
        %       all observed variables other hyper parameters, and 
        %       returns two state structsobjects. 
        %
        %       Sd:     The dynamic state that captures the variable 
        %               values that would be updated at each iteration.
        %               In a sampling program, samples will be collected
        %               based on this state.
        %
        %       Sc:     The constant state that captures fixed quantities,
        %               which might include observation, hyper parameters,
        %               and pre-computed quantities.
        %
        %       S0:     A struct/object that gives (part of) the initial 
        %               states.
        %
        %       optype: The kind of operation to do.
        %               - 'sample':     MCMC / Gibbs sampling
        %               - 'varinfer':   Variational inference        
        %               This argument is supplied by the driver function.
        %
        
        [Sd, Sc] = update(prg, Sd, Sc);
        % Updates the states
        %
        %   [Sd, Sc] = prg.update(Sd, Sc);
        %
        %       This method updates the states. From the view of a driver
        %       function, in each iteration, the update method will get
        %       invoked once.
        %
        %       In typical cases, the update might only modify Sd, and 
        %       keep Sc fixed. However, Sc is allowed to be updated in
        %       proper cases, e.g. the program might make modification
        %       of Sc for performance optimization, etc.
        %
        
        Sp = make_output(prg, Sd, Sc);
        % Makes the output sample from states
        %
        %   Sp = prg.make_output(Sd, Sc);
        %
        %       This method extracts the relevant quantities from the
        %       states and makes an output sample.
        %
        %       In many cases, Sp can be a struct comprised a subset
        %       of fields from Sd. For example, in a mixture model 
        %       program, Sc might comprise inferred/estimated variable
        %       values along with log-likelihood tables, etc, and 
        %       only those variable values are output to Sp.
        %
        %       It is advisable to make Sp in form of a column vector,
        %       or a struct scalar, as in several important driver
        %       functions, collected samples will be concatenated
        %       horizontally into a sequence.
        %
        
        objv = evaluate_objective(prg, Sd, Sc);
        % Evaluates the objective value of the current state
        %
        %   objv = prog.evaluate_objective(Sd, Sc);
        %       
        %       Evaluates the objective value based on current states.
        %       
        %
        
    end
    
end

