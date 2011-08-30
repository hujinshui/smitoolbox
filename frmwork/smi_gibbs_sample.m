function R = smi_gibbs_sample(program, opts, initcfgs)
% Performs (generalized) Gibbs sampling over a smi framework
%
%   R = smi_gibbs_sample(program, opts);
%   R = smi_gibbs_sample(program, opts, initcfgs);
%
%       Performs Gibbs sampling over a statistical modeling and inference
%       framework. 
%
%       Input arguments:
%       - program:      An SMI program upon which the sampling is 
%                       performed.
%
%       - opts:         The options that control the sampling procedure.
%                       This can be made by calling gs_options.
%
%       - initcfgs:     The initial configuration(sample) to begin with.
%
%                       If only one MCMC chain is run, then initcfgs must
%                       be a struct with each field corresponding to a
%                       variable in the framework.
%
%                       If one wants to run multiple MCMC chains in 
%                       parallel, initcfgs can be a struct array, with 
%                       each element providing an initial configuration.
%
%                       If proper initialization steps have been added
%                       to the program, then initcfgs can be omitted.
%
%       Output arguments:
%       - R:            The cell array comprised of collected samples. 
%                       Each cell of R is a sequence of samples collected 
%                       from an MCMC chain, in form of a struct array.
%                       Each element of such an array is a sample, whose
%                       fields correspond to output variables.
%
%                       The size of R is nr x n0, where nr is the number
%                       of repeating experiments for each initial config,
%                       which is specified by the option opts.nrepeat;
%                       n0 is the numel(initcfgs).
%
%   Remarks
%   --------
%       1. If there are multiple init-configs, or there are multiple runs
%          for each config, then parfor is used to perform each run in a
%          parallel manner. In this case, one can config a pool of
%          computation sessions properly to make the scheduling efficient.
%

% Created by Dahua Lin, on Aug 25, 2011
%

%% Verify input arguments

if ~isa(program, 'smi_program')
    error('smi_gibbs_sample:invalidarg', ...
        'program should be an object of class smi_program.');
end

bmap = program.block_map;
if ~isfield(bmap, 'default')
    error('smi_gibbs_sample:invalidarg', ...
        'default block is lacking in the program.');
end

blk_default = program.blocks(bmap.default);

blk_init = [];
if isfield(bmap, 'init')
    blk_init = program.blocks(bmap.init);
end

opts = gs_options(opts);
nr = opts.nrepeats;

if nargin < 3
    initcfgs = [];
else
    if ~isstruct(initcfgs)
        error('smi_gibbs_sample:invalidarg', 'initcfgs should be a struct.');
    end        
end
n0 = numel(initcfgs);


%% Main 

if n0 <= 1    
        
    if nr == 1  % single config, single run each
        R = do_sim(program, blk_default, blk_init, opts, initcfgs, 0);
        R = { R };
        
    else        % single config, multi run each
        R = cell(nr, 1);
        
        parfor i = 1 : nr
            R{i} = do_sim(program, blk_default, blk_init, opts, initcfgs, i);
        end
    end
    
else            
    if nr == 1  % multiple configs, single run each
                
        R = cell(1, n0);
        
        parfor i = 1 : n0
            R{i} = do_sim(program, blk_default, blk_init, opts, initcfgs(i), i);
        end        
        
    else        % multiple configs, multi run each
        
        Cfgs = cell(nr, n0);
        for j = 1 : n0
            for i = 1 : nr
                Cfgs{i, j} = initcfgs(j);
            end
        end
        
        R = cell(nr, n0);                
        
        N = nr * n0;
        parfor i = 1 : N
            R{i} = do_sim(program, blk_default, blk_init, opts, Cfgs{i}, i);
        end
    end
    
end


%% core simulation procedure

function S = do_sim(program, blk_default, blk_init, opts, cfg0, ithread)
% single run with a single config

% establish context

if ~opts.check
    ctx = smi_exec_context(program);
else
    ctx = smi_exec_context(program, 'check');
end

displevel = opts.display;

% initialization

if displevel >= 2
    fprintf('[[thread %d]] initializing ...\n', ithread);
end

if ~isempty(cfg0)
    ctx.set_vars(cfg0);
end

if ~isempty(blk_init)
    exec_block(ctx, blk_init);
end

% run sampling

nsamples = opts.nsamples;
S = cell(nsamples, 1);

if displevel < 4
    
    % burn in
    
    if displevel >= 2
        fprintf('[[thread %d]] burning in ...\n', ithread);
    end    
    
    exec_block(ctx, blk_default, opts.burnin);
    
    % main iterations
    
    if displevel >= 2
        fprintf('[[thread %d]] collecting samples ...\n', ithread);
    end
    
    for i = 1 : nsamples
        exec_block(ctx, blk_default, opts.cps);
        S{i} = ctx.extract_vars(opts.output_vars);
        
        if displevel >= 3
            fprintf('[[thread %d]]   %d/%d samples collected.\n', ...
                ithread, i, nsamples);
        end
    end
    
else % displevel >= 4 (iters)
    
    % burn in
    
    fprintf('[[thread %d]] burning in ...\n', ithread);
    
    for t = 1 : opts.burnin
        fprintf('[[thread %d]]     burn-in iter %d/%d\n', ...
            ithread, t, opts.burnin);
        
        exec_block(ctx, blk_default);       
    end
    
    % main iterations
    
    fprintf('[[thread %d]] collecting samples ...\n', ithread);
    
    cps = opts.cps;
    for i = 1 : nsamples
        for t = 1 : cps
            fprintf('[[thread %d]]     sample[%d] iter %d/%d\n', ...
                ithread, i, t, cps);
            
            exec_block(ctx, blk_default);
        end
        S{i} = ctx.extract_vars(opts.output_vars);
        
        fprintf('[[thread %d]]   %d/%d samples collected.\n', ...
            ithread, i, nsamples);        
    end
    
end

% combine collected samples

S = vertcat(S{:});

if displevel >= 1
    fprintf('[[thread %d]] finished.\n', ithread);
end



    