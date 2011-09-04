function R = gmm_demo(optype)
% A script to demo the inference over Gaussian mixture model
%
%   gmm_demo sample;
%   gmm_demo varinfer;
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 14, 2010
%       - Modified by Dahua Lin, on Aug 32, 2011
%       - 

%% verify input

optype = lower(optype);
if ~(strcmp(optype, 'sample') || strcmp(optype, 'varinfer'))
    error('gaussgm_demo:invalidarg', 'Invalid option.');
end    


%% model configuration

d = 2;
n = 1000;
K = 3;

pri_mu = [6; 3];
pri_sig = 8;

gpri0 = gaussd.from_mp(pri_mu, pdmat('s', d, pri_sig^2));

Cx = pdmat([1, 1]');

%% generate data

U0 = gpri0.sample(K);

Gx0 = gaussd.from_mp(U0, Cx);
X = Gx0.sample(n * ones(1, K), 1:K);

%% build program

Cu = gpri0.C;
gprg = gmm_std(d, Cx, Cu);

obs.X = X;
obs.K = K;
obs.mu = 0;

%% run inference

switch optype
    
    case 'sample'

        nsamples = 20;
        
        opts = mcmc_options([], ...
            'burnin', 100, ...
            'nsamples', nsamples, ...
            'ips', 15, ...
            'display', 'sample');
        
        R = smi_mcmc(gprg, obs, [], opts);
        R = R{1};
        
        Us = cat(3, R.U);  % => d x K x nsamples
        Us = permute(Us, [1 3 2]);  % => d x nsamples x K
        
        Zs = vertcat(R.Z);  % => nsamples x N
        Zm = mode(Zs, 1);
        
    case 'varinfer'
        
        opts = varinfer_options([], ...
            'maxiters', 100, ...
            'display', 'eval');
        
        R = smi_varinfer(gprg, obs, [], opts);     
        Us = R.sol.U;
        [~, Zm] = max(R.sol.Q, [], 1);
end

visualize_results(K, X, Us, Zm);


%% Sub functions


function visualize_results(K, X, Us, Zm)

gm = intgroup(K, Zm);

hfig = figure;
set(hfig, 'Name', 'GMM Demo');

set(hfig, 'Position', [0 0 1024, 512]);
movegui(hfig, 'center');

subplot('Position', [0.05, 0.1, 0.43, 0.8]);
plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 5);
hold on;
if ndims(Us) == 2
    plot(Us(1,:), Us(2,:), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
else
    plot(Us(1,:), Us(2,:), 'r+', 'MarkerSize', 12);
end
axis equal;


subplot('Position', [0.52, 0.1, 0.43, 0.8]);

colors = {'r', 'g', 'b', 'm', 'c', 'k'};
for k = 1 : K   
    cr = colors{mod(k-1, numel(colors)) + 1};
    hold on;
    plot(X(1, gm{k}), X(2, gm{k}), [cr '.'], 'MarkerSize', 8);
end

axis equal;

