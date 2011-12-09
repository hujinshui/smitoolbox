function demo_gmrf_sample(imsiz, w, n, intv)
% Demos the use of Gaussian MRF sampling 
%
%   demo_gmrf_sample(imsiz, w, n, intv);
%
%       demonstrates the sampling of images from a Gaussian MRF.
%
%       Input:
%       - imsiz:        the image size, in form of [height, width]
%       - w:            the weights of the spatial links
%       - n:            the number of samples to generate
%       - intv:         the interval between two samples to be collected.
%

% Created by Dahua Lin, on Dec 9, 2011
%

%% main

% model construction

disp('Constructing the MRF model ...');

imh = imsiz(1);
imw = imsiz(2);

r = 1;      % the neighborhood range
gr = gr_local([imh imw], r);
W = gr_wmat(gr, w);

G = agmrf(W, 0.1);  % obtain the G-MRF


% make the sampler

bdim = 10;

y0 = 1:bdim:imh;
x0 = 1:bdim:imw;

y1 = [y0(2:end) - 1, imh];
x1 = [x0(2:end) - 1, imw];

d = imh * imw;
inds = reshape(1:d, imh, imw);

blocks = cell(numel(y0), numel(x0));
for j = 1 : numel(x0)
    for i = 1 : numel(y0)
        v = inds(y0(i):y1(i), x0(j):x1(j));
        blocks{i, j} = v(:);
    end
end
nblocks = numel(blocks);

sampler = gmrf_blk_gibbs(G, blocks);

% simulation

disp('Simulating the chain ...');

s = repmat(1:nblocks, 1, intv);

figure;

x = randn(d, 1);

himg = imshow(get_vimage(x, imh, imw));
title('t = 0');

disp('Press any key to continue move ..');
pause;

for i = 1 : n    
    x = sampler.update(x, s);    
    
    vimg = get_vimage(x, imh, imw);
    set(himg, 'CData', vimg);
    title(sprintf('t = %d', i * intv));
    
    disp('Press any key to continue move ..');
    pause;
end


%% sub functions

function im = get_vimage(x, imh, imw)

I = reshape(x, imh, imw);
I = I + 0.5;
I(I < 0) = 0;
I(I > 1) = 1;
I = im2uint8(I);
im = cat(3, I, I, I);

