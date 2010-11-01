function S = test_gmrf_lbp()
% a simple script to test the implementation of LBP on G-MRF
%

% construct the graph

n = 3;

J11 = eye(2) * 3;
J22 = eye(3) * 4;
J33 = eye(4) * 5;
J12 = rand(2, 3);
J13 = rand(2, 4);

gm = gaussmrf({J11, J22, J33}, [1 2; 1 3], {J12, J13});
dims = gm.dims;

% ground-truth

J0 = full(gm.information_matrix);
C0 = inv(J0);
h0 = rand(sum(dims), 1);
u0 = C0 * h0; %#ok<MINV>

Cs0 = cell(n, 1);
us0 = cell(n, 1);
hs = cell(n, 1);

b = 0;
for i = 1 : n    
    si = b + (1:dims(i));
    
    Cs0{i} = C0(si, si);
    us0{i} = u0(si);
    hs{i} = h0(si);
    
    b = b + dims(i);
end


% run bp

tbp = gmrf_tbp(gm, 1);
tbp.infer_Js();

tbp.initialize_hs(hs);
tbp.infer_hs();

% extract results

Cs = cell(n, 1);
us = cell(n, 1);

for i = 1 : n
    Cs{i} = inv(tbp.r_Js{i});
    us{i} = Cs{i} * tbp.r_hs{i};
end

c_diffs = zeros(1, n);
u_diffs = zeros(1, n);
for i = 1 : n
    c_diffs(i) = norm(Cs{i} - Cs0{i});
    u_diffs(i) = norm(us{i} - us0{i});
end

% output

S.gm = gm;
S.tbp = tbp;
S.Cs0 = Cs0;
S.Cs = Cs;
S.us0 = us0;
S.us = us;
S.c_diffs = c_diffs;
S.u_diffs = u_diffs;

