function R = smi_mcmc(prg, obs, S0, opts)
% Drive a MCMC sampling procedure based on an SMI program
%
%   R = smi_mcmc(prg, obs, S0, opts);
%
%       Drives a MCMC sampling procedure and collects samples.
%
%       Input arguments:
%       - prg:      the SMI program object, which should be an object
%                   of a class derived from smi_prg.
%
%       - obs:      The observations to be passed to the program.
%
%       - S0:       The initial state to be passed to the program.
%
%       - opts:     the MCMC control options, which can be obtained by
%                   calling mcmc_options.
%                   (See the help of mcmc_options for details).
%       
%       Output arguments:
%       - R:        a cell array of sample sequences. 
%
%                   In general, each Markov chain being run will result
%                   in a sequence of samples. Each sample sequence here
%                   is formed by horizontally concatenating the results
%                   yielded by the make_output method of prg.
%
%                   R is a cell array of size nchains x 1, where nchains
%                   if the number of independent chains spawned by the
%                   function, and R{i} is the sample sequence obtained
%                   from the i-th chain.
%
%   Remarks
%   -------
%       - If opts.nrepeats > 0, multiple independent Markov chains
%         will be spawned, running in an parallel way (using parfor).
%         A proper setting of matlabpool would help to enhance the
%         parallel efficiency.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 4, 2011
%

%% verify input arguments

if ~isa(prg, 'smi_prg')
    error('smi_mcmc:invalidarg', 'prg should be a SMI program.');
end

opts = mcmc_options(opts);

%% main

nr = opts.nrepeats;

if nr == 1
    R = {[]};
    R{1} = do_mcmc(prg, obs, S0, opts, 0);
else
    R = cell(1, nr);
    parfor i = 1 : nr
        R{i} = do_mcmc(prg, obs, S0, opts, i);
    end
end


%% core MCMC function

function sseq = do_mcmc(prg, obs, S0, opts, ithread)

displevel = opts.display;

% initialization

if displevel >= 2
    fprintf('[[thread %d]] initializing ...\n', ithread);
end

[Sd, Sc] = prg.initialize(obs, S0, 'sample');
    
% burn in

if displevel >= 2
    fprintf('[[thread %d]] burning in ...\n', ithread);
end

for t = 1 : opts.burnin
    if displevel >= 4
        fprintf('[[thread %d]]     burn-in iter %d/%d\n', ...
            ithread, t, opts.burnin);
    end    
    [Sd, Sc] = prg.update(Sd, Sc);       
end

% main iterations

if displevel >= 2
    fprintf('[[thread %d]] collecting samples ...\n', ithread);
end

nsamples = opts.nsamples;
sseq = cell(1, nsamples);

ips = opts.ips;
for i = 1 : nsamples
    for t = 1 : ips
        if displevel >= 4
            fprintf('[[thread %d]]     sample[%d] iter %d/%d\n', ...
                ithread, i, t, ips);
        end        
        [Sd, Sc] = prg.update(Sd, Sc);
    end
    sseq{i} = prg.make_output(Sd, Sc);
    
    if displevel >= 3
        fprintf('[[thread %d]]   %d/%d samples collected.\n', ...
            ithread, i, nsamples);
    end
end
    
% combine samples into a sequence

sseq = horzcat(sseq{:});

if displevel >= 1
    fprintf('[[thread %d]] finished.\n', ithread);
end


