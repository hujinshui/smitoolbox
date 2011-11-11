function gsm_demo(method)
% A simple demo of Gaussian scale mixture fitting
%
%   gsm_demo(method);
%
%   Here method can be either 'ip' or 'em'.
%

% Created by Dahua Lin, on Nov 2, 2011
%

%% configure model

K = 4;
sigma = [1, 4, 12, 36].';
p = [0.4, 0.2, 0.25, 0.15].';

ext = 100;

%% generate data

n = 1e6;
z = ddsample(p, n);

grps = intgroup(K, z);

x = zeros(1, n);
for k = 1 : K
    g = grps{k};
    x(g) = randn(1, numel(g)) * sigma(k);
end

xi = linspace(-ext, ext, 5000);

%% model fitting

tol = 1e-2;
tolx = 1e-4;

[p_e, s_e] = gsm_fit(x, [], K, 'method', method, 'tol', tol, 'tolx', tolx, ...
    'display', 'iter');

[s_e, sc] = sort(s_e, 1, 'ascend');
p_e = p_e(sc);

display(p_e);
display(s_e);


%% plotting

pdf_gt = gsm_pdf(xi, p, sigma);
pdf_es = gsm_pdf(xi, p_e, s_e);

figure;
semilogy(xi, pdf_gt, xi, pdf_es);
legend({'ground-truth', 'estimated'});
xlabel('x');
ylabel('pdf (in log-scale)');


function P = gsm_pdf(x, p, sigma)

E = exp( ( -0.5 ./ (sigma.^2) ) * (x.^2) );    
P = (p ./ sqrt((2 * pi) * (sigma.^2)))' * E;


