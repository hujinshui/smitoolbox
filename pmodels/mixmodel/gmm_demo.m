function [R, Gs] = gmm_demo(optype, cf, op)
% A script to demo the inference over Gaussian mixture model
%
%   R = gmm_demo('sample',   cf);
%   R = gmm_demo('varinfer', cf);
%
%       Runs a Gaussian mixture model demo using the specified inference
%       methodology. The covariance matrix has the form specified by
%       cf, which can be either 's' (scale form), 'd' (diagonal form), 
%       or 'f' (full matrix form).
%
%       The output is the struct that captures inference result.
%
%   R = gmm_demo('sample',   cf, 'tied_cov');
%   R = gmm_demo('varinfer', cf, 'tied_cov');
%
%       Runs the demo under the setting that the covariance matrix
%       is shared across all mixture components.
%
%   [R, Gs] = gmm_demo( ...);
%
%       Additionally returns the resultant Gaussian model (of class
%       gaussd).
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 14, 2010
%       - Modified by Dahua Lin, on Aug 31, 2011
%       - Modified by Dahua Lin, on Sep 28, 2011
%

%% verify input

optype = lower(optype);
if ~(strcmp(optype, 'sample') || strcmp(optype, 'varinfer'))
    error('gmm_demo:invalidarg', 'Invalid option.');
end    

if ~(ischar(cf) && isscalar(cf) && (cf == 's' || cf == 'd' || cf == 'f'))
    error('gmm_demo:invalidarg', 'Invalid cf.');
end

if nargin < 3
    c_tied = false;
else
    if ~(ischar(op) && strcmp(op, 'tied_cov'))
        error('gmm_demo:invalidarg', 'The 3rd argument is invalid.');
    end
    c_tied = true;
end


%% model configuration

d = 2;
n = 1000;
K = 3;

pri_mu = [6; 3];
pri_sig = 4;

gpri0 = gaussd.from_mp(pri_mu, pdmat('s', d, pri_sig^2), 'ip');


switch cf
    case 's'
        alpha = 5;
        beta = 5;
        cpri0 = invgammad(1, alpha, beta, 'pre');
        Vx = cpri0.sample(K);
        Cx = pdmat(cf, d, Vx);

    case 'd'
        alpha = 5;
        beta = 5;
        cpri0 = invgammad(d, alpha, beta, 'pre');
        Vx = cpri0.sample(K);
        Cx = pdmat(cf, d, Vx);
        
    case 'f'      
        Phi = pdmat('s', d, 15);
        wdeg = 15;
        cpri0 = invwishartd(Phi, wdeg, 'pre');
        Cx = cpri0.sample(K);        
        
end

%% generate data

U0 = gpri0.sample(K);

Gx0 = gaussd.from_mp(U0, Cx);
X = Gx0.sample(n * ones(1, K), 1:K);

%% build program

if ~c_tied
    gm = gaussgm(d);
else
    gm = gaussgm(d, 'tied_cov');
end

% gprg = fmm_std(gm, {gpri0, cf}, K);
gprg = fmm_std(gm, {gpri0, cpri0}, K);
gprg.dalpha = 3;

%% run inference

switch optype
    
    case 'sample'

        nsamples = 30;
        
        opts = mcmc_options([], ...
            'burnin', 100, ...
            'nsamples', nsamples, ...
            'ips', 15, ...
            'display', 'sample');
        
        R = smi_mcmc(gprg, X, [], opts);
        R = R{1};
        
        all_params = [R.params];  % => d x K x nsamples
        Us = cat(3, all_params.U);
        Um = mean(Us, 3);
        Cs = [all_params.Cx];
        if cf == 's' || cf == 'd'
            Vs = cat(3, Cs.v);
            Vm = mean(Vs, 3);
            Cm = pdmat(cf, d, Vm);
        else
            Vs = cat(4, Cs.v);
            Vm = mean(Vs, 4);
            Cm = pdmat('f', d, Vm);
        end
        
        Gs = gaussd.from_mp(Um, Cm);
        
        Zs = vertcat(R.Z);  % => nsamples x N
        Zm = mode(Zs, 1);
        
    case 'varinfer'
        
        opts = varinfer_options([], ...
            'maxiters', 500, ...
            'display', 'eval');
        
        R = smi_varinfer(gprg, X, [], opts);     
        params = R.sol.params;
        [~, Zm] = max(R.sol.Q, [], 1);
        Gs = gaussd.from_mp(params.U, params.Cx);
        
end

visualize_results(K, X, Gs, Zm);


%% Sub functions


function visualize_results(K, X, Gs, Zm)

gm = intgroup(K, Zm);

hfig = figure;
set(hfig, 'Name', 'GMM Demo');

set(hfig, 'Position', [0 0 1024, 512]);
movegui(hfig, 'center');

subplot('Position', [0.05, 0.1, 0.43, 0.8]);
plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 5);
hold on;
Gs.plot_ellipse(2, 'r-', 'LineWidth', 2);
axis equal;

subplot('Position', [0.52, 0.1, 0.43, 0.8]);

colors = {'r', 'g', 'b', 'm', 'c', 'k'};
for k = 1 : K   
    cr = colors{mod(k-1, numel(colors)) + 1};
    hold on;
    plot(X(1, gm{k}), X(2, gm{k}), [cr '.'], 'MarkerSize', 8);
end

axis equal;

