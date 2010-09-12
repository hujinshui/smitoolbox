% Compute safe dot product between vectors in A and B
%
%   v = safedot(A, B);
%       computes the "safe" dot product between vectors in A amd B.
%       The safe dot product between vectors x and y is defined as follows
%   
%           sum_i x(i) * y(i)
%
%       Here, when either x(i) or y(i) is zero, then x(i) * y(i) is 
%       regarded as zero, no matter whether the other one is inf or nan.
%
%       By default, the computation is along the first non-singleton
%       dimension.
%
%   v = safedot(A, B, dim);
%       computes the "safe" dot product between vectors in A and B along
%       specified dimension.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 12, 2010
%

