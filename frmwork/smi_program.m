classdef smi_program
    % The class to capture an statistical modeling and inference program
    %
    %   Each program is a collection of compiled code blocks, which 
    %   can refer to each other. A program can be run by an smi driver.
    %
    %   Each block is a struct with the following fields:
    %      
    %   - name:     the name of the block
    %   - length:   the (first-level) length of execution sequence (n)
    %   - ops:      an n x 1 char vector of operation types
    %   - ids:      an n x 1 vector of integer reference ids
    %   - nrs:      an n x 1 vector of repeating-time numbers
    %   - ivs:      an n x 1 cell array of input variable ids
    %   - ovs:      an n x 1 cell array of output variable ids
    %
    %   Each step (say the i-th one) is characterized by 
    %   - the operation type (ops(i)): 
    %       ops(i) == 'f':  a function execution
    %       ops(i) == 'b':  a reference to other block
    %   - the reference id (ids(i)):
    %       if ops(i) == 'f': this is the id of the function to call
    %       if ops(i) == 'b': this is the id of the block to call
    %   - the input variable ids (ivs{i}): (applies only when op == 'f')
    %   - the output variable ids(ovs{i}): (applies only when op == 'f')
    %
        
    % Created by Dahua Lin, on Aug 29, 2011
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        vars;           % the table of variables
        funcs;          % the table of functions
        blocks;         % the table of blocks
        
        var_map;        % variable name -> id
        func_map;       % function name -> id
        block_map;      % block name -> id
    end
    
    properties(Dependent)
        num_vars;
        num_funcs;
        num_blocks;
    end
    
    methods
        function n = get.num_vars(prg)
            n = numel(prg.vars);
        end
        
        function n = get.num_funcs(prg)
            n = numel(prg.funcs);
        end
        
        function n = get.num_blocks(prg)
            n = numel(prg.blocks);
        end
    end    
    
    %% Compilation
    
    methods(Static)
        
        function prg = compile(frmwork)
            % Compiles an SMI framework into a program
            %
            %   prg = smi_program.compile(frmwork);
            %
            %       This function constructs an SMI program by compiling
            %       an SMI framework.
            %
            
            % compile vars
            
            vs = frmwork.vars;
            vmap = smi_program.make_map(vs);
            
            fs = frmwork.funcs;
            fmap = smi_program.make_map(fs);            
            
            % compile blocks
            
            B = frmwork.blocks;            
            nb = frmwork.num_blocks;            
            blk_map = smi_program.make_map(B);
                        
            blks = repmat( struct( ...
                'name', [], 'length', [], ...
                'ops', [], 'ids', [], 'nrs', [], 'ivs', [], 'ovs', [] ...
                ), nb, 1);
            
            for i = 1 : nb
                blks(i) = compile_block(B{i}, vmap, fmap, blkmap);
            end
            
            prg = smi_program;
            
            prg.vars = vs;
            prg.funcs = fs;
            prg.blocks = blks;
            
            prg.var_map = vmap;
            prg.func_map = fmap;
            prg.block_map = blk_map;
        end
        
    end
    
    
    methods(Static, Access='private')
        
        function cblk = compile_block(blk, funcs, varmap, funmap, blkmap)
            % Compiles one block
            
            cblk.name = blk.name;
            
            n = numel(blk.steps);
            cblk.length = n;
            
            cblk.ops = char(n, 1);
            cblk.ids = zeros(n, 1);
            cblk.nrs = zeros(n, 1);
            cblk.ivs = cell(n, 1);
            cblk.ovs = cell(n, 1);
            
            for i = 1 : n                
                if strcmp(s.type, 'f')
                    
                    cblk.ops(i) = 'f';
                    
                    f_id = smi_program.get_entity_id( ...
                        funmap, s.func_name, 'function name');
                    cblk.ids(i) = f_id;
                    
                    cblk.nrs(i) = s.nr;
                    
                    fo = funcs(f_id);
                    cblk.ivs{i} = smi_program.get_var_ids( ...
                        fo.input_slot_names, s.input_slots, s.input_vars, varmap);
                    cblk.ovs{i} = smi_program.get_var_ids( ...
                        fo.output_slot_names, s.output_slots, s.output_vars, varmap);
                
                elseif strcmp(s.type, 'b')
                    
                    cblk.ops(i) = 'b';
                    
                    b_id = smi_program.get_entity_id( ...
                        blkmap, s.block_name, 'block name');
                    cblk.ids(i) = b_id;
                    
                    cblk.nrs(i) = s.nr;
                    cblk.ivs{i} = [];
                    cblk.ovs{i} = [];
                                        
                else
                    error('smi_program:compile_error', ...
                        'Invalid step type %s', s.type);
                end
            end                                    
        end
        
        
        function M = make_map(L)
            % Makes name->id map from a list of named items
            
            M = [];
            n = numel(L);
            for i = 1 : n
                M.(L{i}.name) = i;
            end            
        end
        
        
        function id = get_entity_id(M, name, title)
            % Gets the id of an entity from a map
            
            if isfield(M, name)
                id = M.(name);
            else
                error('There is no %s %s in the program.', title, name);
            end
        end
        
        function vids = get_var_ids(funname, slotnames, slots, vars, vmap)
            % Extract the ids to form an array of vars
                        
            vids = zeros(1, numel(slotnames));
            for i = 1 : numel(svpairs)
                s = slots{i};
                v = vars{i};
                
                [tf, islot] = ismember(s, slotnames);
                if ~tf
                    error('There is no slot in function %s called %s', ...
                        funname, s);
                end
                
                [tf, ivar] = ismember(v, vmap);
                if ~tf
                    error('There is no variable name %s in the program', v);
                end
               
                if vids(islot) ~= 0
                    error('The slot %s of function %s is connected more than once.', ...
                        s, funname);
                end
                vids(islot) = ivar;
            end
        end
        
        
    end
    
    
end
