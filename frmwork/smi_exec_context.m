classdef smi_exec_context < handle
    % The Execution context for running SMI program
    %
    
    % Created by Dahua Lin, on Aug 30, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')
        program;                % the program being executed
        values;                 % the current values of all variables
        is_checked = false;     % whether to check variable validities
    end
    
    %% Construction and Manipulation
    
    methods
        
        function ctx = smi_exec_context(program, op)
            % Establishes the context for running an SMI program
            %
            %   ctx = smi_exec_context(program);
            %
            %       establishes the context for running the input program.
            %
            %   ctx = smi_exec_context(program, 'check');
            %
            %       establishes a context which verifies the inputs and
            %       outputs at each step.
            %
            
            ctx.program = program;
            
            nv = program.num_vars;
            ctx.values = cell(nv, 1);
            
            if nargin >= 2
                if ~(ischar(op) && strcmpi(op, 'check'))
                    error('smi_exec_context:invalidarg', ...
                        'The 2nd argument can only be ''check''.');
                end
                ctx.is_checked = true;
            end
            
        end
        
        
        function set_vars(ctx, V)
            % Sets a set of values to variables in the context
            %
            %   ctx.set_vars(V);
            %
            %       Here, V is a struct, such that field names
            %       are variable names, and the value of each field
            %       is to be set to the corresponding variables
            %       in the context.
            %
            
            if ~(isstruct(V) && isscalar(V))
                error('smi_exec_context:invalidarg', ...
                    'The input to set_vars should be a struct scalar.');
            end
            
            vmap = ctx.program.var_map;
                        
            vns = fieldnames(V);
            n = numel(vns);
            for i = 1 : n
                vn = vns{i};
                if isfield(vmap, vn)
                    vid = vmap.(vn);
                    ctx.values{vid} = V.(vn);
                else
                    error('smi_exec_context:invalidarg', ...
                        'The variable %s is not in the program context.', vn);
                end
            end            
        end
        
        
        function V = extract_vars(ctx, vnames)
            % Extracts a set of variables from context to form a struct
            %
            %   V = extract_vars(ctx);
            %   V = extract_vars(ctx, vnames);
            %       
            %       It extracts the values of the variables (whose names
            %       are given by vnames) to form a struct V.
            %
            %       If vnames is empty or omitted, it extracts all
            %       variables. 
            %
            
            if nargin < 2 || isempty(vnames)
                vnames = {ctx.program.vars.name};
            else
                if ~iscellstr(vnames)
                    error('smi_exec_context:invalidarg', ...
                        'vnames should be either empty or a cell array of strings.');
                end
            end
            
            vmap = ctx.program.var_map;
            
            n = numel(vnames);
            for i = 1 : n
                vn = vnames{i};
                if isfield(vmap, vn)
                    vid = vmap.(vn);
                    V.(vn) = ctx.values{vid};
                else
                    error('smi_exec_context:invalidarg', ...
                        'The variable %s is not in the program context.', vn);
                end
            end            
        end        
        
    end
    
    
    %% Execution
    
    methods
    
        function exec_step(ctx, s)
            % Execute a step in the context
            %
            %   ctx.exec_step(s);
            %
            %       Executes the step s in context. The variables
            %       captured by the context are updated from the 
            %       execution outputs.
            %
            
            prg = ctx.program;
            chk = ctx.is_checked;            
            
            if s.op(1) == 'f'
                
                nr = s.nr;
                for i = 1 : nr
                    f = prg.funcs(s.id);
                    fobj = f.obj;
                    
                    in_vals = cell(1, f.n_in);
                    out_vals = cell(1, f.n_out);
                    
                    % prepare input                    
                    for j = 1 : f.n_in
                        iv = s.ivs(j);
                        if iv > 0
                            in_vals{j} = ctx.values{iv};
                        end
                    end
                    
                    % check the validity of inputs
                    if chk
                        ctx.check_vars(true, s, f, in_vals)
                    end
                    
                    % evaluate function                    
                    [out_vals{:}] = fobj.evaluate(s.ovs > 0, in_vals{:});
                    
                    % check the validity of outputs                    
                    if chk
                        ctx.check_vars(false, s, f, out_vals)
                    end
                                        
                    % dispatch output                    
                    for j = 1 : f.n_out
                        ov = s.ovs(j);
                        if ov > 0
                            ctx.values{ov} = out_vals{j};
                        end
                    end
                    
                end
                
            elseif s.op(1) == 'b'
                
                nr = s.nr;
                for i = 1 : nr
                    exec_block(ctx, prg.blocks(s.id));
                end
                
            else
                error('Invalid step opchar %s', s.type);
            end
            
        end
        
        
        function exec_block(ctx, b, nr)
            % Execuate a block in the context
            %
            %   ctx.exec_block(b);
            %   ctx.exec_block(b, nr);
            %
            %       Executes the block b in the context. nr is the
            %       number of repeats. If omitted, nr is set to 1.
            %
            
            if nargin < 3
                nr = 1;
            end
            
            if nr <= 0
                return;
            end
            
            ns = numel(b.steps);
            if nr == 1
                for i = 1 : ns
                    exec_step(ctx, b.steps(i));
                end
            else
                for t = 1 : nr
                    for i = 1 : ns
                        exec_step(ctx, b.steps(i));
                    end
                end
            end            
        end
            
    end
    
    
    methods(Access='private')
        
        function check_vars(ctx, is_input, step, f, vals)
            % check the validity of input/ouput variables
                                    
            vs = ctx.program.vars;
            
            if is_input
                
                in_slots = f.obj.input_slots;
                n_in = numel(in_slots);
                if numel(vals) ~= n_in
                    error('smi_exec_context:rterror', ...
                        'The set size of input vars is incorrect.');
                end
                
                ivs = step.ivs;
                
                for j = 1 : n_in
                    if ~ivs(j); continue; end
                    vv = vals{j};
                    if ~smi_verify_var(vv, in_slots(j))
                        error('smi_exec_context:rterror', ...
                            'Invalid value: var %s => func %s.(%s).', ...
                            vs(ivs(j)).name, in_slots(j).name, f.name);
                    end
                end
                
            else % is output
                
                out_slots = f.obj.output_slots;
                n_out = numel(out_slots);
                if numel(vals) ~= n_out
                    error('smi_exec_context:rterror', ...
                        'The set size of output vars is incorrect.');
                end
                
                ovs = step.ovs;
                
                for j = 1 : n_out
                    if ~ovs(j); continue; end
                    vv = vals{j};
                    if ~smi_verify_var(vv, out_slots(j))
                        error('smi_exec_context:rterror', ...
                            'Invalid value: func %s.(%s) => var %s.', ...
                            out_slots(j).name, f.name, vs(ovs(j)).name);
                    end
                end
                
            end            
        end
        
        
    end
    
    
end


