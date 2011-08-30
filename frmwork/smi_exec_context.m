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
    
    %% Construction and 
    
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
            
            if s.type(1) == 'f'
                
                nr = s.nr;
                for i = 1 : nr
                    f = prg.funcs(s.id);
                    
                    in_vars = cell(1, f.n_in);
                    out_vars = cell(1, f.n_out);
                    
                    % prepare input
                    
                    for j = 1 : f.n_in
                        iv = s.ivs(j);
                        if iv > 0
                            in_vars{j} = ctx.values{iv};
                        end
                    end
                    
                    % evaluate function
                    
                    [out_vars{:}] = f.obj.evaluate(f.ovs > 0, in_vars{:});
                    
                    % dispatch output
                    
                    for j = 1 : f.n_out
                        ov = s.ovs(j);
                        if ov > 0
                            ctx.values{ov} = out_vars{j};
                        end
                    end
                    
                end                
                
            elseif s.type(1) == 'b'
                
                nr = s.nr;                
                for i = 1 : nr
                    exec_block(ctx, prg.blocks(s.id));
                end
                
            else
                error('Invalid step type %s', s.type);
            end
            
        end
        
        
        function exec_block(ctx, b)
            % Execuate a block in the context
            %
            %   ctx.exec_block(b);
            %
            %       Executes the block b in the context.
            %
            
            ns = numel(b.steps);
            for i = 1 : ns
                exec_step(ctx, b.steps(i));
            end
            
        end
        
    
    end
    
    
end
