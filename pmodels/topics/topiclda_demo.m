function topiclda_demo()
%TOPICLDA_DEMO A simple program to demo the use of LDA
%
%   R = TOPICLDA_DEMO;
%

%% Configuration

V = 50;     % the number of words in the vocabulary
K = 3;      % the number of topics
n0 = 500;   % the number of training documents
n1 = 1000;  % the number of testing documents

nw = 1000;  % the number of words per document


%% Generate data

disp('Generating document data ...');

[alpha0, Beta0] = gen_cfg(V, K);

C0 = gen_corpus(alpha0, Beta0, nw, n0);
C1 = gen_corpus(alpha0, Beta0, nw, n1);


%% Train model

disp('Training LDA model ...');

eta = 0.01;
S = topiclda_em(V, eta);
S = S.initialize(C0, [], K);

opts = varinfer_options([], 'maxiters', 200, 'tol', 1e-6, 'display', 'off');
S = varinfer_drive(S, opts);

Beta = S.sol.Beta;
alpha = S.sol.alpha;

% permute topic labels

D = zeros(K, K);
for k = 1 : K
    for l = 1 : K
        D(k, l) = ddkldiv(Beta0(:,k), Beta(:,l));
    end
end

[~, tmap] = min(D, [], 1);
Beta(:, tmap) = Beta;
D(:, tmap) = D;
alpha(tmap) = alpha;

%% Show results

disp('============================');

disp('Topic prior:');
disp('  ground-truth:  '); disp(alpha0.' / sum(alpha0));
disp('  solved-result: '); disp(alpha.' / sum(alpha)); 

disp('Divergence between True Topics and Solved Topics');
disp(D);

figure;
for k = 1 : K
    subplot(K,1,k);
    bar([Beta0(:,k), Beta(:,k)]);
end

disp('Evaluating Perplexity values ...');

ppx0_tr = calc_ppx(alpha0, Beta0, C0);
ppx0_te = calc_ppx(alpha0, Beta0, C1);
ppx1_tr = calc_ppx(alpha, Beta, C0);
ppx1_te = calc_ppx(alpha, Beta, C1);

[alpha_ram, Beta_ram] = gen_cfg(V, K);
ppxr_tr = calc_ppx(alpha_ram, Beta_ram, C0);
ppxr_te = calc_ppx(alpha_ram, Beta_ram, C1);


disp('                        True      Estimated        Random');
fprintf('  On training:  %12.4f   %12.4f   %12.4f\n', ppx0_tr, ppx1_tr, ppxr_tr);
fprintf('  On testing:   %12.4f   %12.4f   %12.4f\n', ppx0_te, ppx1_te, ppxr_te);


disp(' ');


%% Data generation function

function [alpha, Beta] = gen_cfg(V, K)

alpha = 1 + rand(K, 1) * 2;
Beta = exp(randn(V, K) * 1.5) * 2;
Beta = bsxfun(@times, Beta, 1 ./ sum(Beta, 1));


function C = gen_corpus(alpha, Beta, nw, n)

[V, K] = size(Beta);
T = dird_sample(alpha, n);   

% draw topics
Ts = ddsample(T, nw);

% draw words

gs = intgroup(K, Ts(:));
words = zeros(nw, n);

for k = 1 : K
    cg = gs{k};
    ws = ddsample(Beta(:,k), numel(cg));
    words(cg) = ws;
end

% count words

C = zeros(V, n);
for i = 1 : n
    C(:,i) = intcount(V, words(:,i));
end


function v = calc_ppx(alpha, Beta, C)

K = size(Beta, 2);
Ginit = bsxfun(@plus, alpha, sum(C, 1) / K);
Gam = topiclda_varinfer(Beta, alpha, C, [], Ginit, 500, 1e-8);
v = topiclda_ppx(Beta, C, Gam);


