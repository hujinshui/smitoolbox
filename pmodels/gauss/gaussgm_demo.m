function R = gaussgm_demo()
% A script to demo the inference over Gaussian generative model
%
%   gaussgm_demo;
%

% Created by Dahua Lin, on Nov 14, 2010
% Modified by Dahua Lin, on Aug 32, 2011


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

% gpri = gaussd.from_mp(pri_mu, pdmat('s', d, 0.01^2));
gpri = gpri0;
gprg = gmm_std(X, K, Cx, gpri, [], 'gs');

%% run inference

nsamples = 20;

opts = gs_options([], ...
    'burnin', 100, ...
    'nsamples', nsamples, ...
    'cps', 15, ...
    'output_vars', {'U', 'Z'}, ...
    'check', true, ...
    'display', 'sample');

R = smi_gibbs_sample(gprg, opts);
R = R{1};

%% extract results

Us = cat(3, R.U);  % => d x K x nsamples
Us = permute(Us, [1 3 2]);  % => d x nsamples x K 

Zs = vertcat(R.Z);  % => nsamples x N 
Zm = mode(Zs, 1);
gm = intgroup(K, Zm);

%% visualization

hfig1 = figure;
set(hfig1, 'Name', 'GMM Demo (Inferred Models)');

plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 5);

for k = 1 : K
    hold on;
    plot(Us(1,:,k), Us(2,:,k), 'r+', 'MarkerSize', 12);    
end
axis equal;


hfig2 = figure;
set(hfig2, 'Name', 'GMM Demo (Inferred Labels)');

colors = {'r', 'g', 'b', 'm', 'c', 'k'};
for k = 1 : K   
    cr = colors{mod(k-1, numel(colors)) + 1};
    hold on;
    plot(X(1, gm{k}), X(2, gm{k}), [cr '.'], 'MarkerSize', 8);
end

axis equal;


