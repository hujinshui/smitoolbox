function ker = ggkernel2d(type, hsiz, param)
% Generate a kernel for gridgraph2d
%
%   ker = ggkernel2d('unif', hsiz);
%       
%       generates a uniform kernel of specified size. The siz argument
%       can be in form of [m, n], or just a scalar n, which is eqivalent
%       to [n n]. It represents the half-size. In the output, the size of
%       ker is [m+1, n+1], representing a kernel of actual size (2m+1) x
%       (2n+1).
%
%   ker = ggkernel2d('ellip', hsiz);
%
%       generates a uniform kernel whose masked region is of ellipse shape.
%
%   ker = ggkernel2d('gauss', hsiz, sigma);
%
%       generates a Gaussian kernel. The argument sigma can be a scalar 
%       for isotropic case, or it can be a pair for using different scales
%       along different dimensions.
%
%   Remarks
%   -------
%       - the generated kernel has ker(1, 1) = 0. 
%

% Created by Dahua Lin, on Feb 8, 2011
%

%% verify input

if ~ischar(type)
    error('ggkernel2d:invalidarg', 'type should be a char string.');
end

if ~(isnumeric(hsiz) && numel(hsiz) <= 2)
    error('ggkernel2d:invalidarg', 'hsiz should be a numeric scalar or pair');
end

if isscalar(hsiz)
    m = hsiz;
    n = hsiz;
else
    m = hsiz(1);
    n = hsiz(2);
end

%% main

switch lower(type)
    case 'unif'
        if nargin >= 3
            error('ggkernel2d:invalidarg', 'unif kernel does not have a parameter.');
        end
        ker = constmat(m+1, n+1, 1);
        ker(1) = 0;
        
    case 'ellip'
        if nargin >= 3
            error('ggkernel2d:invalidarg', 'ellip kernel does not have a parameter.');
        end
        [xx, yy] = meshgrid((0:n)/n, (0:m)/m);
        rr = xx.^2 + yy.^2;
        ker = double(rr <= 1);
        ker(1) = 0;
        
    case 'gauss'
        sigma = param;
        if ~(isfloat(sigma) && numel(sigma) <= 2)
            error('ggkernel2d:invalidarg', ...
                'sigma should be either a numeric scalar or a pair.');
        end
        if isscalar(sigma)
            sx = sigma(1);
            sy = sigma(1);
        else
            sx = sigma(2);
            sy = sigma(1);
        end
        [xx, yy] = meshgrid((0:n) / sx, (0:m) / sy);
        ker = exp(-0.5 * (xx.^2 + yy.^2));
        ker(1) = 0;
        
    otherwise
        error('ggkernel2d:invalidarg', ...
            'Invalid kernel type %s', type);
end




