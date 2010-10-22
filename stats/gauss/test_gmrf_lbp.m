function S = test_gmrf_lbp()
% a simple script to test the implementation of LBP on G-MRF
%

% construct the graph

n = 8;
edges = [1 2; 2 3; 2 4; 3 5; 3 6; 4 7; 4 8];
db = 3;

G = construct_random_graph(n, db, edges, 0.5);

% do inference

bp = blk_gmrf_bp(G);

max_iter = 100;
for t = 1 : max_iter
    bp.update_msg(1:n);
    bp.update_msg(n:-1:1);
end

% get results

J = full(G.precision_matrix);
h = G.potential_vector;

[mu_s, sigma_s] = get_result(bp);

mu = vertcat(mu_s{:});
var_s = cellfun(@diag, sigma_s, 'UniformOutput', false);
sigma = vertcat(var_s{:});

% output

S.G = G;
S.bp = bp;
S.J = J;
S.h = h;
S.mu0 = J \ h;
S.sigma0 = diag(inv(J));

S.mu = mu;
S.sigma = sigma;


function G = construct_random_graph(n, db, edges, a)

d = n * db;

W = zeros(d, d);
ne = size(edges, 1);

for i = 1 : n
    Wi = rand(db, db);
    Wi = 0.5 * (Wi + Wi.');
    W((i-1)*db + (1:db), (i-1)*db+(1:db)) = Wi;
end

for i = 1 : ne
    
    s = edges(i, 1);
    t = edges(i, 2);
    
    Wi = rand(db, db);
    Wi = 0.5 * (Wi + Wi.');
    
    W((s-1)*db + (1:db), (t-1) * db + (1:db)) = Wi;
    W((t-1)*db + (1:db), (s-1) * db + (1:db)) = Wi';    
end

ds = sum(W, 1);
J = diag(ds + a) - W;
h = randn(d, 1) * 3;

Js = cell(n, n);
hs = cell(n, 1);

for i = 1 : n
    hs{i} = h((i-1)*db+(1:db));
    Js{i,i} = J((i-1)*db+(1:db), (i-1)*db+(1:db));
end

for i = 1 : ne
    
    s = edges(i, 1);
    t = edges(i, 2);
   
    Js{s,t} = J((s-1)*db+(1:db), (t-1)*db+(1:db));
    Js{t,s} = J((t-1)*db+(1:db), (s-1)*db+(1:db));
end

G = blk_gmrf(hs, Js);

assert(isequal(full(G.precision_matrix), J));
assert(isequal(G.potential_vector, h));


