classdef gmrf_blk_gibbs
    % The class that implements a block gibbs sampler for Gaussian MRF
    %
    
    % Created by Dahua Lin, on Dec 9, 2011
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        gmrf;       % the underlying gaussian model
             
        blocks;     % the cell array of index vectors for all blocks
        
        Jmats;      % the cell array of partial precision matrices
        Lmats;      % the cell array of Cholesky matrices  
        Fmats;      % the cell array of forward matrices
    end
    
    
    %% constructor
    
    methods
        
        function obj = gmrf_blk_gibbs(G, blocks)
            % Construct a block gibbs sampler for a Gaussian MRF
            %
            %   obj = gmrf_blk_gibbs(G, blocks);
            %
            %   Input parameters:
            %   G:      the underlying Gaussian MRF. G should be a 
            %           gaussd object, with G.ty == 'c'.
            %
            %   blocks: a cell array of blocks, each cell is a vector
            %           of indices of variables within a block.            
            %
            
            % verify arguments
            
            if ~(is_gaussd(G) && G.ty == 'c' && G.n == 1)
                error('gmrf_blk_gibbs:invalidarg', ...
                    'G should be a gaussd struct with G.ty == ''c'' and G.n == 1.');
            end
            
            if ~iscell(blocks)
                error('gmrf_blk_gibbs:invalidarg', ...
                    'blocks should be a cell array of index vectors.');
            end
            
            nb = numel(blocks);
            for i = 1 : nb
                b = blocks{i};
                if ~(isvector(b) && isnumeric(b))
                    error('gmrf_blk_gibbs:invalidarg', ...
                        'blocks{%d} is invalid.', i);
                end
            end            
            
            % block-wise data
            
            Js = cell(nb, 1);
            Ls = cell(nb, 1);
            Fs = cell(nb, 1);
            J = G.J.v;
            for i = 1 : nb
                b = blocks{i};
                Ji = full(J(b, b));
                Js{i} = Ji;
                Ls{i} = chol(Ji);
                
                Fi = J(b, :);
                Fi(:, b) = [];
                Fs{i} = Fi;
            end            
            
            % set fields
            
            obj.gmrf = G;
            obj.blocks = blocks;
            
            obj.Jmats = Js;
            obj.Lmats = Ls;      
            obj.Fmats = Fs;
        end
                
    end
    
    
    %% sampling
    
    
    methods
       
        function X = update(obj, X, s)
            % Performs Gibbs updates 
            %
            %   X = obj.update(X, s);
            %
            %       Updates given samples following a sequence of steps.                   
            %       Each step re-draws the values of a specific block.
            %
            %       Inputs:
            %       - X:        the samples to be updated. It can be
            %                   a column vector of size d x 1, or a 
            %                   matrix of size d x n (n is the number
            %                   of chains that are being simulated).
            %
            %       - s:        the sequence of blocks to be updated.
            %                   The same block can repeat in the sequence.
            %
            
            % take useful fields
            
            d = obj.gmrf.d;
            h = obj.gmrf.h;
            B = obj.blocks;
            Js = obj.Jmats;
            Ls = obj.Lmats;
            Fs = obj.Fmats;
            
            lopts.UT = true;
            
            % verify inputs
            
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == d)
                error('gmrf_blk_gibbs:invalidarg', ...
                    'X should be a real matrix with size(X,1) == gmrf.d.');
            end
            n = size(X, 2);
            
            if ~(isvector(s) && isnumeric(s) && isreal(s))
                error('gmrf_blk_gibbs:invalidarg', ...
                    's should be sequence of block indices.');
            end
            if size(s, 1) > 1; s = s.'; end
            
            % Go!
            
            for i = s
                
                b = B{i};
                di = numel(b);
                             
                Ji = Js{i};
                Li = Ls{i};     
                Fi = Fs{i};
                                
                Xr = X;
                Xr(b, :) = [];
                dH = Fi * Xr;
                
                % Hi = h(b) - Fi * Xr;
                if isequal(h, 0)
                    Hi = -dH;
                else
                    if n == 1
                        Hi = h(b) - dH;
                    else
                        Hi = bsxfun(@minus, h(b), -dH);
                    end
                end
                Ui = Ji \ Hi;
                               
                Z = randn(di, n);
                dX = linsolve(Li, Z, lopts);
                
                X(b, :) = Ui + dX;
            end
        end
        
    end
    
    
end

