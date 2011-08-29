function R = smi_gibbs_sample(frmwork, opts, initcfgs)
% Performs (generalized) Gibbs sampling over a smi framework
%
%   R = smi_gibbs_sample(frmwork, opts, initcfgs);
%
%       Performs Gibbs sampling over a statistical modeling and inference
%       framework. 
%
%       Input arguments:
%       - frmwork:      The framework upon which sampling is performed.
%                       It should be an object of the class smi_frmwork.
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

if ~isa(frmwork, 'smi_frmwork')
    error('smi_gibbs_sample:invalidarg', ...
        'frmwork should be an object of class smi_frmwork.');
end

opts = gs_options(opts);
nr = opts.nrepeats;

if ~isstruct(initcfgs)
    error('smi_gibbs_sample:invalidarg', 'initcfgs should be a struct.');
end
n0 = numel(initcfgs);

if n0 == 1    
    repr = preprocess_cfg(frmwork, initcfgs);
    
    if nr == 1  % single config, single run each
        R = do_sim(frmwork, opts, repr, 0);
        R = { R };
        
    else        % single config, multi run each
        R = cell(nr, 1);
        
        parfor i = 1 : nr
            R{i} = do_sim(frmwork, opts, repr, i);
        end
    end
    
else            
    if nr == 1  % multiple configs, single run each
        
        reprs = cell(1, n0);
        for j = 1 : n0
            reprs{j} = preprocess_cfg(frmwork, initcfgs(j));
        end
        
        R = cell(1, n0);
        
        parfor i = 1 : n0
            R{i} = do_sim(frmwork, opts, reprs{i}, i);
        end        
        
    else        % multiple configs, multi run each
        
        reprs = cell(nr, n0);
        for j = 1 : n0
            repr = preprocess_cfg(frmwork, initcfgs(j));            
            for i = 1 : nr
                reprs{i, j} = repr;
            end            
        end
        
        R = cell(nr, n0);                
        
        N = nr * n0;
        parfor i = 1 : N
            R{i} = do_sim(frmwork, opts, reprs{i}, i);
        end
    end
    
end


%% core simulation procedure

function S = do_sim(frmwork, opts, repr, ithread)
% single run with a single config

% extract fields

seq0 = frmwork.init_seq;
len0 = numel(seq0);

cseq = frmwork.cycle_seq;
clen = numel(cseq);

fs = frmwork.funcs;
insts0 = frmwork.instructions0;
insts = frmwork.instructions;
vmap = frmwork.var_map;

T0 = opts.burnin;
N = opts.nsamples;
nc = opts.cps;
displevel = opts.display;

if displevel >= 1
    if ithread > 0
        dprefix = ['[[thread ' int2str(ithread), ']] '];
    else
        dprefix = '';
    end
end

% get ids of variable to be output to samples
if isempty(opts.output_vars)
    ovars = {frmwork.vars.name};
else
    ovars = opts.output_vars;
end
ovar_ids = zeros(size(ovars));
for i = 1 : numel(ovars)
    ovar_ids(i) = vmap.(ovars{i});
end

% initialize result container
S = cell(1, N);

if displevel >= 2
    disp([dprefix, 'simulation started']);
end

% initialization 

if len0 > 0
    for i = 1 : len0

        % pick the current step
        st = insts0(i);
        f = fs(st.f_id).func;

        if displevel >= 4
            disp([dprefix, sprintf('    step %d: %s', ...
                i, fs(st.f_id).name)]);
        end

        % form inputs
        vin = smi_from_repr(repr, st.v_in);

        % do update computation
        vout = cell(1, st.n_out);
        [vout{:}] = f.evaluate(st.out_flags, vin{:});

        % put the results back to repr
        for j = 1 : st.n_out
            if st.v_out(j) > 0
                repr{st.v_out(j)} = vout{j};
            end
        end
    end

    if displevel >= 2
        disp([dprefix, 'initialization finished']);
    end
end


% main loop
Tstart = -T0;
Tend = (nc - 1) * N;

c = 0; % sample counter            
for t = Tstart : Tend

    if t == 0 && displevel >= 2
        disp([dprefix, 'burn-in finished']);
    end

    if displevel >= 3                    
        disp([dprefix, sprintf('  cycle %d ...', t)]);
    end

    % do update cycle
    for i = 1 : clen 

        % pick the current step
        st = insts(i);
        f = fs(st.f_id).func;                  

        if displevel >= 4
            disp([dprefix, sprintf('    step %d: %s', ...
                i, fs(st.f_id).name)]);
        end

        % form inputs                                        
        vin = smi_from_repr(repr, st.v_in);

        % do update computation
        vout = cell(1, st.n_out);
        [vout{:}] = f.evaluate(st.out_flags, vin{:});

        % put the results back to repr
        for j = 1 : st.n_out
            if st.v_out(j) > 0
                repr{st.v_out(j)} = vout{j};
            end
        end
    end       

    % collect sample                
    if t >= 0 && mod(t, nc) == 0
        c = c + 1;                    
        for i = 1 : numel(ovars)
            sp.(ovars{i}) = repr{ovar_ids(i)};
        end
        S{c} = sp;
    end
end

% convert results to struct form
S = vertcat(S{:});  

if displevel >= 1
    disp([dprefix, 'simulation done']);
end



%% auxiliary functions

function repr = preprocess_cfg(frmwork, cfg)
% form internal representation and verify
        
nv = frmwork.num_vars;
vs = frmwork.vars;

% form internal representation
repr = cell(1, nv);
for i = 1 : nv
    vn = vs(i).name;
    if isfield(cfg, vn)
        repr{i} = cfg.(vn);
    end
end
    

