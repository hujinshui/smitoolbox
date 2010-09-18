% Compute safe dot product between vectors in A and B
%
%   v = safedot(A, B);
%       computes the "safe" dot product between vectors in A amd B.
%       The safe dot product between vectors x and y is defined as follows
%   
%           sum_i x(i) * y(i)
%
%       Here, when either x(i) or y(i) is zero, then x(i) * y(i) is 
%       regarded as zero, even the other one is inf or nan.
%
%       If A and B are both vectors, then it returns the safe dot product
%       between A and B. If either A or B is a row/column vector, then
%       it returns the dot product between that vector and the row/column
%       vectors of the other matrix. If both are non-vector matrices, then
%       it returns the dot product values between the columns in A and B.
%
%   v = safedot(A, B, dim);
%       computes the "safe" dot product between vectors in A and B along
%       specified dimension.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 12, 2010
%

