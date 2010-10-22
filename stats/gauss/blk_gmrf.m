classdef blk_gmrf
    % The class to represent a block-based Gaussian MRF model
    %
    
    % Created by Dahua Lin, on Oct 21, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')        
        nblocks;    % the number of blocks (K)
        dims;       % the dimension of each block [K x 1]
        tdim;       % the total dimension
        
        hs;         % the cell array of potential vectors [K x 1 cells]
        Js;         % the cell array of precision matrices [K x K cells]        
    end
    
    methods
        
        function obj = blk_gmrf(hs, Js)
            % Constructs a block-based Gaussian MRF object
            %
            %   obj = blk_gmrf(hs, Js);
            %       constructs a G-MRF model using information-form
            %       parameterization. 
            %
            %       Suppose there are nb blocks, then hs should be
            %       an nb x 1 cell array, and Js should be an nb x nb
            %       cell array. Such that hs{k} is the potential vector
            %       of the k-th block, and Js{k, k'} is the information
            %       matrix part between k-th and k'-th block.
            %
            %       Note that Js{k, k'} can be empty if there is no
            %       direct edges between the k-th and k'-th blocks.
            %
            
            if ~(iscell(hs) && isvector(hs))
                error('blk_gmrf:invalidarg', 'hs should be a cell vector.');
            end
            if size(hs, 2) > 1
                hs = hs.';
            end
            
            if ~(iscell(Js) && ndims(Js) == 2)
                error('blk_gmrf:invalidarg', 'Js should be a cell matrix.');
            end
            
            K = numel(hs);
            if ~isequal(size(Js), [K K])
                error('blk_gmrf:invalidarg', ...
                    'The sizes of hs and Js are inconsistent.');
            end
            
            ds = zeros(K, 1);
            for k = 1 : K
                h_k = hs{k};
                if ~(isfloat(h_k) && ndims(h_k) && size(h_k, 2) == 1)
                    error('blk_gmrf:invalidarg', ...
                        'each h_k should be a numeric column vector.');
                end
                ds(k) = size(h_k, 1);
            end
            
            for k = 1 : K
                for l = 1 : K
                    dk = ds(k);
                    dl = ds(l);
                    
                    J_kl = Js{k, l};
                    if ~isempty(J_kl)
                        if ~(isfloat(J_kl) && isequal(size(J_kl), [dk dl]))
                            error('blk_gmrf:invalidarg', ...
                                'Some part of information matrix is invalid.');
                        end
                    end                    
                end
            end
            
            obj.nblocks = K;
            obj.dims = ds;
            obj.tdim = sum(ds);
            obj.hs = hs;
            obj.Js = Js;            
        end
        
        
        function h = potential_vector(obj)
            % Get the full potential vector of the model
            %
            %   h = obj.potential_vector;
            %
            
            h = vertcat(obj.hs{:});            
        end
        
        
        function J = precision_matrix(obj)
            % Get the full precision matrix of the model
            %
            %   J = obj.precision_matrix;
            %
            % Note: the resultant matrix is a sparse matrix.
            %
            
            K = obj.nblocks;
            bs = [0; cumsum(obj.dims(1:end-1))];
        
            i_s = cell(K, K);
            j_s = cell(K, K);
            w_s = cell(K, K);
            
            for k = 1 : K
                for l = 1 : K                    
                    J_kl = obj.Js{k, l};
                    if ~isempty(J_kl)                    
                        [i,j,w] = find(J_kl);
                        
                        i_s{k, l} = i + bs(k);
                        j_s{k, l} = j + bs(l);
                        w_s{k, l} = w;
                    end                    
                end
            end
            
            ii = vertcat(i_s{:});
            jj = vertcat(j_s{:});
            ww = vertcat(w_s{:});
            
            td = obj.tdim;            
            J = sparse(ii, jj, ww, td, td);            
        end
        
    end

end



