function R = mcmc_drive(S, opts)
%MCMC_DRIVE Drive the running of an MCMC procedure.
%
%   R = MCMC_DRIVE(S, opts);
%
%       Drives a MCMC sampling procedure and collects samples.
%
%       Input arguments:
%       - S0:       The SMI state (an instance of class smi_state) that
%                   has been initialized.
%
%       - opts:     The MCMC control options, which can be obtained by
%                   calling mcmc_options.
%                   (See the help of mcmc_options for details).
%       
%       Output arguments:
%       - R:        A cell array of samples (each sample
%                   results from calling the output method of the state).
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

if ~isa(S, 'smi_state')
    error('mcmc_drive:invalidarg', 'S should be an smi_state object.');
end

if ~S.is_ready()
    error('varinfer_drive:invalidarg', 'S has not been ready.');
end

opts = mcmc_options(opts);

%% main

displevel = opts.display;

% burn in

if displevel >= 2
    fprintf('Burning in ...\n');
end

for t = 1 : opts.burnin
    if displevel >= 4
        fprintf('     burn-in iter %d/%d\n', t, opts.burnin);
    end    
    S = update(S);
end

% main iterations

if displevel >= 2
    fprintf('Collecting samples ...\n');
end

nsamples = opts.nsamples;
R = cell(1, nsamples);

ips = opts.ips;
for i = 1 : nsamples
    if ips == 1
        S = update(S);
    else
        for t = 1 : ips
            if displevel >= 4
                fprintf('     sample[%d] iter %d/%d\n', i, t, ips);
            end        
            S = update(S);
        end
    end
    R{i} = output(S);
    
    if displevel >= 3
        fprintf('   %d/%d samples collected.\n', i, nsamples);
    end
end
    


