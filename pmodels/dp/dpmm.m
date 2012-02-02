classdef dpmm < smi_prg
    % The class to implement a DP mixture model (DPMM) program
    %
    
    % Created by Dahua Lin, on Sep 17, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')   
        amodel;     % the underlying atom model (of class genmodel_base)
        basedist;   % the base distribution
        alpha;      % the concentration parameter    
    end
    
    % configuration
    
    properties
        dead_tolratio = 0.3;    % the tolerable ratio of dead atoms
    end    
        
    
    %% Constructor
    
    methods
        
        function prg = dpmm(model, basedist, alpha)
            % Creates an DPMM program
            %
            %   sol = dpmm(model, basedist, alpha, ...);
            %
            %       Creates a DPMM program, given the underlying
            %       atom model (model), the base distribution, and 
            %       the concentration parameter (alpha).
            %
            
            % verify inputs
            
            if ~isa(model, 'genmodel_base')
                error('dpmm:invalidarg', ...
                    'The input model should be of a class derived from genmodel_base.');
            end
            
            if ~model.is_supported_optype('atom')
                error('dpmm:invalidarg', ...
                    'The input model does not support the optype ''atom''.');
            end
            
            if ~model.is_valid_prior(basedist)
                error('dpmm:invalidarg', ...
                    'basedist is not a valid prior w.r.t to the atom model.');
            end
            
            if ~(isfloat(alpha) && isscalar(alpha) && isreal(alpha) && alpha > 0)
                error('dpmm:invalidarg', ...
                    'alpha should be a positive real number.');
            end                        
                        
            % set fields
            
            prg.amodel = model;
            prg.basedist = basedist;
            prg.alpha = alpha;
        end
        
    end
       
    
    %% Interface methods
    
    methods
                        
        function [sol, Sc] = initialize(prg, obs, S0, optype)
            % Initialize the states given observations and relevant info
            %
            %   [sol, Sc] = prg.initialize(obs, S0, 'sample');
            %
            %       The method initialize the solution and a struct that
            %       captures relevant information about the observation
            %       and inherited atoms.
            %
            %       S0 is either empty (use default settings), or 
            %       a struct with all or part of the following fields
            %       that customize the settings.
            %
            %       - inherits:     an inherit struct that captures the
            %                       atoms inherited from prior source.
            %                       (can be made by dpmm_inherits)
            %                       If omitted, there is no inheritance.
            %
            %       - initcap:      the initial capacity of the solution.
            %
            %       - sol:          the initial solution (of class dpmm_sol).
            %
            %       In the output, 
            %       - sol:      the initial solutio (of class dpmm_sol)
            %       - Sc:       a struct of static information:
            %                   - obs:          the observation array
            %                   - inherits:     the inheritance struct
            %                   - nh;           the number of inherits
            %
            
            % basics
            
            amdl = prg.amodel;
            n = amdl.get_num_observations(obs);
            
            H = [];
            c0 = [];
            sol = [];            
            
            % check init fields
            
            if ~isempty(S0)                
                if ~(isstruct(S0) && isscalar(S0))
                    error('dpmm:invalidarg', ...
                        'S0 should be either empty or a struct scalar.');
                end
                
                if isfield(S0, 'inherits') && ~isempty(S0.inherits)
                    H = S0.inherits;
                    if ~(isstruct(H) && isfield(H, 'tag') && ...
                            strcmp(H.tag, 'dpmm_inherits'))
                        error('dpmm:invalidarg', ...
                            'The inherits field is not valid.');
                    end                                        
                end
                
                if isfield(S0, 'initcap') && ~isempty(S0.initcap)
                    c0 = S0.initcap;
                    if ~(isnumeric(c0) && isscalar(c0) && c0 >= 0)
                        error('dpmm:invalidarg', ...
                            'initcap should be a non-negative scalar.');
                    end
                end
                
                if isfield(S0, 'sol') && ~isempty(S0.sol)
                    sol = S0.sol;
                    if ~isa(sol, 'dpmm_sol')
                        error('dpmm:invalidarg', ...
                            'sol should be an instance of class dpmm_sol.');
                    end
                    if sol.nobs ~= n
                        error('dpmm:invalidarg', ...
                            '#obs is inconsistent with sol.nobs.');
                    end
                end                
            end
            
            if ~(ischar(optype) && strcmpi(optype, 'sample'))
                error('dpmm:invalidarg', ...
                    'The only supported operation type is ''sample''.');
            end
            
            % construct solution
            
            if isempty(sol)
                if isempty(c0)
                    sol = dpmm_sol(amdl, prg.basedist, obs);
                else
                    sol = dpmm_sol(amdl, prg.basedist, obs, c0);
                end
            end
            
            % construct Sc
            
            Sc.inherits = H;
            Sc.obs = obs;
            if isempty(H)
                Sc.nh = 0;
            else
                Sc.nh = H.num;
            end            
        end
        
        
        function [sol, Sc] = update(prg, sol, Sc)
            % Updates the solution
            %
            %   [sol, Sc] = prg.update(sol, Sc);
            %
            %       Updates the solution by running re-sampling steps.
            %
            
            amdl = prg.amodel;            
            dist0 = prg.basedist;
            obs = Sc.obs;
            
            sol = sol.update_labels(amdl, dist0, obs, prg.alpha);
            sol = sol.update_atoms(amdl, dist0, obs, Sc.inherits);       
            sol = sol.prune(Sc.nh, prg.dead_tolratio);
        end
        
        
        function S = make_output(prg, sol, ~) %#ok<MANU>
            % Makes the output sample from states
            %
            %   S = make_output(prg, sol, Sc);
            %
            %       This returns a struct that representa s sample,
            %       which contains the following fields:
            %
            %       - natoms;
            %       - max_atom_id;
            %       - atom_ids;
            %       - atoms;
            %       - atom_counts;
            %       - iatoms;
            %       - labels := atom_ids(iatoms)
            % 
            
            S.natoms = sol.natoms;
            S.max_atom_id = sol.max_atom_id;
            S.atom_ids = sol.atom_ids;
            S.atoms = sol.atoms;
            S.atom_counts = sol.atom_counts;
            S.iatoms = sol.iatoms;
            S.labels = sol.labels;            
        end
        
        
        function objv = evaluate_objective(prg, sol, Sc) %#ok<STOUT,MANU,INUSD>
            % Evaluate objective of a solution (not implemented)
            %
                        
            error('dpmm:notimplemented', ...
                'Evaluation of objective is not implemented for DPMM.');
        end
        
    end
       
    
end


