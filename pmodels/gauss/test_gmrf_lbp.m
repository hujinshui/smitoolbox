function S = test_gmrf_lbp()
% a simple script to test the implementation of LBP on G-MRF
%

% construct the graph

J11 = eye(2) * 3;
J22 = eye(3) * 4;
J33 = eye(4) * 5;
J12 = rand(2, 3);
J13 = rand(2, 4);

gm = gaussmrf({J11, J22, J33}, [1 2; 1 3], {J12, J13});


% set bp

tbp = gmrf_tbp(gm, 1);



% output

S.gm = gm;
S.tbp = tbp;