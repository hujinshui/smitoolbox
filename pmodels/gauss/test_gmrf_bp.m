function S = test_gmrf_bp(name)
% a simple script to test the implementation of LBP on G-MRF
%
%   

% construct the model

gm = make_gmrf(name);
h0 = rand(sum(gm.dims), 1);
[Cs0, us0, hs] = compute_ground_truth(gm, h0);

% main 

tbp = run_tbp(gm, hs);
tbp_r = process_results(tbp, Cs0, us0);

lbp = run_lbp(gm, hs);
lbp_r = process_results(lbp, Cs0, us0);

% output

S.gm = gm;
S.Cs0 = Cs0;
S.us0 = us0;
S.tbp_r = tbp_r;
S.lbp_r = lbp_r;


%% sub functions

function gm = make_gmrf(name)

switch name
    case 'toy1'
        J11 = eye(2) * 3;
        J22 = eye(3) * 4;
        J33 = eye(4) * 5;
        J12 = rand(2, 3);
        J13 = rand(2, 4);

        gm = gaussmrf({J11, J22, J33}, [1 2; 1 3], {J12, J13});
        
    case 'toy2'
        Js = cell(8, 1);
        Js{1} = randpsd([4, 3]);
        Js{2} = randpsd([6, 3, 4]);
        Js{3} = randpsd([5, 2.5, 3]);
        Js{4} = randpsd([3, 3]);
        Js{5} = 4;
        Js{6} = randpsd([3, 4.2, 3.5]);
        Js{7} = randpsd([4, 3]);
        Js{8} = 5;
        
        J12 = rand(2, 3);
        J23 = rand(3, 3);
        J34 = rand(3, 2);
        J35 = rand(3, 1);        
        J36 = rand(3, 3);
        J27 = rand(3, 2);
        J18 = rand(2, 1);
        
        gm = gaussmrf(Js, ...
            [1 2; 2 3; 3 4; 3 5; 3 6; 2 7; 1 8], ...
            {J12, J23, J34, J35, J36, J27, J18});
        
        J = full(gm.information_matrix);
        min_e = min(eig(J));
        
        if min_e < 1;
            for i = 1 : numel(Js)
                di = size(Js{i}, 1);
                Js{i} = Js{i} + eye(di) * (1 - min_e);
            end
        end
        
    otherwise
        error('Unknown model name %s', name);        
end


function [Cs0, us0, hs] = compute_ground_truth(gm, h0)

n = gm.nnodes;
J0 = full(gm.information_matrix);
C0 = inv(J0);
u0 = C0 * h0; %#ok<MINV>

Cs0 = cell(n, 1);
us0 = cell(n, 1);
hs = cell(n, 1);

dims = gm.dims;
b = 0;
for i = 1 : n    
    si = b + (1:dims(i));
    
    Cs0{i} = C0(si, si);
    us0{i} = u0(si);
    hs{i} = h0(si);
    
    b = b + dims(i);
end



function tbp = run_tbp(gm, hs)

tbp = gmrf_tbp(gm, 1);
tbp.infer_Js();

tbp.initialize_hs(hs);
tbp.infer_hs();


function lbp = run_lbp(gm, hs)

lbp = gmrf_lbp(gm);
lbp.infer_Js('tol', 1e-12);

lbp.initialize_hs(hs);
lbp.infer_hs('tol', 1e-12);


function r = process_results(bp, Cs0, us0)

n = numel(Cs0);

Cs = cell(n, 1);
us = cell(n, 1);

for i = 1 : n
    Cs{i} = inv(bp.r_Js{i});
    us{i} = Cs{i} * bp.r_hs{i};
end

c_diffs = zeros(1, n);
u_diffs = zeros(1, n);
for i = 1 : n
    c_diffs(i) = norm(Cs{i} - Cs0{i});
    u_diffs(i) = norm(us{i} - us0{i});
end

r.sol = bp;
r.Cs = Cs;
r.us = us;
r.c_diffs = c_diffs;
r.u_diffs = u_diffs;


function M = randpsd(eigvals)
% generating random positive definite matrix

n = numel(eigvals);
R = orth(rand(n, n));
M = R * diag(eigvals) * R';
M = 0.5 * (M + M');




