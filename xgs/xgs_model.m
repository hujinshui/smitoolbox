classdef xgs_model < handle
    % The class to represent a probabilistic model for 
    % eXtensible Gibbs Sampling (XGS)
    %
    
    % Created by Dahua Lin, on Aug 24, 2011
    %
    
    
    %% properties
    
    % public properties 
    
    properties(GetAccess='public', SetAccess='private')
        
        num_vars = 0;       % the number of variables
        num_steps = 0;      % the number of distinct sampling steps
        cycle_seq = 0;      % the sequence of steps performed per cycle
                
        var_names = {};     % the cell array of all variable names
        steps = [];         % the struct array comprised of steps
                                
        is_compiled = false;    % whether the model has been compiled        
    end
    
    % private properties (implementation-specific)
    
    properties(GetAccess='private', SetAccess='private')
        
        var_map;        % the map: var name -> var id
        step_map;       % the map: step name -> step id   
        
        comp_steps;     % the compiled step information
    end
    
    
    %% methods for model building
    
    methods
        
        function add_var(model, name)
            % Add a named variable to the model
            %
            %   model.add_var(name);
            %
            %       Adds a variable of the input name to the model.
            %
            %   Note: the variable name must be a valid name of a
            %         matlab variable, i.e. isvarname(name) is true.
            %
            
            if ~isvarname(name)
                error('xgs_model:invalidarg', ...
                    'The input name is not a valid variable name.');
            end
            
            if model.is_compiled
                error('xgs_model:invalidarg', ...
                    'Cannot add variable to a compiled model.');
            end
            
            nv = model.num_vars;
            model.var_names{nv+1, 1} = name;
            model.num_vars = nv+1;                   
        end
        
        
        function add_step(model, name, step, conn)
            % Add a named step to the model
            %
            %   model.add_step(name, step, 
            %       {'slot_name1', 'var_name1'; 
            %        'slot_name2', 'var_name2'; ...} );
            %
            %       Adds a step of the input name to the model.
            %       Each step is uniquely identified by its name, which
            %       must be a valid matlab variable name.
            %       
            %       The argument step must be an object of a class that
            %       is a derived class of xgs_step.
            %
            %       Each step is connected with a set of named variables
            %       via its named slots. The 3rd argument is to set up
            %       such connections between slots and variables.
            %
            
            if ~isvarname(name)
                error('xgs_model:invalidarg', ...
                    'The input name is not a valid step name.');
            end
            
            if model.is_compiled
                error('xgs_model:invalidarg', ...
                    'Cannot add step to a compiled model.');
            end
            
            if ~(iscell(conn) && ~isempty(conn) && ...
                    ndims(conn) == 2 && size(conn,2) == 2)
                error('xgs_model:invalidarg', ...
                    'conn should be a cell array with two columns.');
            end
            
            if ~all(cellfun(@isvarname, conn(:)))
                error('xgs_model:invalidarg', ...
                    'Some element in conn is not a valid name.');
            end                        
            
            ns = model.num_steps;
            model.steps(ns+1, 1) = struct( ...
                'name', name, ...
                'obj', step, ...
                'slots', conn(:, 1), ...
                'vars', conn(:, 2));
            model.num_steps = ns+1;            
        end
        
        
        function set_cycle(model, seq)
            % Set the sequence of steps performed in each cycle
            %
            %   model.set_cycle(seq);
            %
            %       Sets the sequence of steps performed in each cycle
            %       of sampling. 
            %
            %       The input argument seq should be a vector of numbers,
            %       each being the index of a step.
            %
            %       The cycle can be set before and after compilation.
            %
            
            if ~(isnumeric(seq) && isvector(seq) && isfull(seq))
                error('xgs_model:invalidarg', ...
                    'The cycle sequence should be a numeric vector.');
            end
            
            if ~isa(seq, 'int32'); seq = double(seq); end            
            if size(seq, 1) > 0; seq = seq.'; end
            
            if any(seq < 1) || any(seq > model.num_steps)
                error('xgs_model:invalidarg', ...
                    'Some numbers in the input seq are out of boundary.');
            end
            
            model.cycle_seq = seq;            
        end
           
    end
    
    
    %% methods for compilation
    
    methods
        
        function compile(model)
            % Compile the model
            %
            %   compile(model) or model.compile();
            %
            %       Compiles the model such that it can be run more
            %       efficiently during sampling.
            %
            %       If the model has alread been compiled, this 
            %       method is a no-action.
            %
            
            if model.is_compiled
                return;
            end
           
            vmap = xgs_model.build_idmap(model.var_names);
            smap = xgs_model.build_idmap([model.steps.name]);            
            csteps = xgs_model.compile_steps(model.steps, vmap); 
            
            model.var_map = vmap;
            model.step_map = smap;
            model.comp_steps = csteps;
        end 
    end
    
    
    methods(Static, Access='private')
        % private implementation of compilation
                
        function idmap = build_idmap(Vs, itemname)
            nv = numel(Vs);
            
            if numel(unique(Vs)) < nv
                error('xgs_model:comperr', ...
                    'Duplicated %s encountered.', itemname);
            end
            
            for i = 1 : nv
                idmap.(Vs{i}) = i;
            end
        end
        
        
        function csteps = compile_steps(S, vmap)
            
            ns = numel(S);
            csteps = repmat( struct( ...
                's_in', [], ...
                'v_in', [], ...
                's_out', [], ...
                'v_out', [], ...
                'out_flags', []), ...
                ns, 1);
            
            for i = 1 : ns                
                o = S(i).obj;                
                ss = S(i).slots;
                vs = S(i).vars;
                
                m_in = o.num_input_slots();
                m_out = o.num_output_slots();
                                                
                s_in = zeros(1, m_in);
                v_in = zeros(1, m_in);
                s_out = zeros(1, m_out);
                v_out = zeros(1, m_out);
                                                    
                % set in/out map
                                    
                for j = 1 : ni
                   si = o.get_slot_info(ss{j});
                   if ~isempty(si.in_id)
                       s_in(si.in_id) = j;
                       v_in(si.in_id) = vmap.(vs{j});
                   end
                   if ~isempty(si.out_id)
                       s_out(si.out_id) = j;
                       v_out(si.out_id) = vmap.(vs{j});
                   end                       
                end                                
                
                % test whether slots are properly connected
                                
                if ~(o.test_slots(s_in > 0, s_out > 0))
                    error('xgs_model:comperr', ...
                        'The slots of step %s are not properly connected.', ...
                        S(i).name);
                end                                
                
                % store info
                csteps(i).s_in = s_in;
                csteps(i).v_in = v_in;
                csteps(i).s_out = s_out;
                csteps(i).v_out = v_out;
                csteps(i).out_flags = s_out > 0;
            end
            
        end
        
    end
    
    
    %% methods for run
    
    methods
        
        function run(model, opts, initsample)
            
        end
        
        function par_run(model, opts, initsamples)
            
            
        end
                
    end
    
    
end
