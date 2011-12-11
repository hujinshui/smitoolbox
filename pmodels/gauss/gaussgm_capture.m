function [dh, dJ] = gaussgm_capture(X, w, A)
% Capture the observations for Gaussian generative model
%
%   The Gaussian generative model, with parameter u, is formulated as
%   
%       x ~ N(A * u, Cx),   or x ~ N(u, Cx) if A is identity.
%
%   Here, let d be the dimension of x, and q be that of u.
%
%   Here, the observation can be captured by the conjugate update as
%
%       dh = sum_i w_i (A' * Cx * x_i);
%       dJ = sum_i w_i (A' * Cx * A);
%
%
%   [dh, dJ] = gaussgm_capture(X);
%   [dh, dJ] = gaussgm_capture(X, w);
%   [dh, dJ] = gaussgm_capture(X, w, A);
%
%       Evaluates the conjugate updates from the given data.
%
%       Inputs:
%       - X:        the sample matrix, size: d x n
%       - w:        the weights of the samples, size: m x n, or a scalar.
%                   If m > 1, then multiple updates are to be evaluated,
%                   each corresponding to a group in w.
%       - A:        the transform matrix. If omitted, it is identity.
%
%       Outputs:
%       - dh:       the updates to the potential vector [q x m]
%       - dJ:       a pdmat struct with dJ.n == m and dJ.d == q.
%

% Created by Dahua Lin, on Dec 10, 2011
%

%% verify inputs




