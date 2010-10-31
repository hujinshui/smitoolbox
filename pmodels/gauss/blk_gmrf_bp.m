classdef blk_gmrf_bp < handle
    % The class to represent a Belief-propagation solution on block GMRF
    %
    
    % Created by Dahua Lin, on Oct 21, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')        
        nblocks;    % the number of blocks: K
        dims;       % the dimension of each block
        nbs;        % the neighbors of each block
        
        hs;     % the cell array of potential vectors [K x 1 cells]
        Js;     % the cell array of precision matrices [K x K cells]
        
        msg_hs;     % the messages (for potential vector(s))
        msg_Js;     % the messages (for precision matrices)
        
        pre_hs;  
        pre_Js;
    end
    
    methods
        
        function obj = blk_gmrf_bp(G)
            % constructs a BP solution on a block Gaussian MRF
            %
            %   obj = blk_gmrf_bp(G);
            %
            
            if ~isa(G, 'blk_gmrf')
                error('blk_gmrf_bp:invalidarg', 'G should be a blk_gmrf object.');
            end
            
            obj.nblocks = G.nblocks;
            obj.dims = G.dims;
            
            obj.hs = G.hs;
            obj.Js = G.Js;
            
            K = G.nblocks;
            
            emap = ~cellfun(@isempty, G.Js);   
            emap(1 + (K+1)*(0:K-1)) = 0;
            obj.nbs = cell(K, 1);
            for k = 1 : K
                obj.nbs{k} = find(emap(k,:));
            end
            
            m_hs = cell(K, K);
            m_Js = cell(K, K);   
            
            p_hs = cell(K, 1);
            p_Js = cell(K, 1);
            
            % initialize
            
            ds = G.dims;            
            for k = 1 : K                
                p_hs{k} = G.hs{k};
                p_Js{k} = G.Js{k, k};
                
                for l = obj.nbs{k}
                    dl = ds(l); 
                    
                    m_hs{k, l} = zeros(dl, 1);
                    m_Js{k, l} = zeros(dl, dl);
                end                                
            end
                        
            obj.pre_hs = p_hs;
            obj.pre_Js = p_Js;
            
            obj.msg_hs = m_hs;
            obj.msg_Js = m_Js;
        end
    
        
        function update_msg(obj, nodes)
            % Updates the message sent from the selected nodes
            %
            %   obj.update_msg(nodes);
            %
            
            hs0 = obj.hs;
            Js0 = obj.Js;
            
            m_hs = obj.msg_hs;
            m_Js = obj.msg_Js;
            
            if size(nodes, 1) > 1
                nodes = nodes.';
            end
            
            for i = nodes            
                Ni = obj.nbs{i};
                
                if isempty(Ni); return; end
                
                p_h = hs0{i} + sum(horzcat(m_hs{Ni, i}), 2);
                p_J = Js0{i,i} + sum(cat(3, m_Js{Ni, i}), 3);
                
                obj.pre_hs{i} = p_h;
                obj.pre_Js{i} = p_J;
                
                for j = Ni
                    Jn = p_J - m_Js{j, i};
                    hn = p_h - m_hs{j, i};
                    
                    J_ij = Js0{i, j};
                    J_ji = Js0{j, i};
                    
                    mm = - (J_ji * (Jn \ [J_ij hn]));
                    
                    m_J = mm(:,1:end-1);
                    m_J = (m_J + m_J') * 0.5;
                    m_h = mm(:, end);
                    
                    obj.msg_hs{i, j} = m_h;
                    obj.msg_Js{i, j} = m_J;
                end            
            end
        end
        
        
        function [mu, sigma] = get_result(obj, i)
            % Extracts the result at current stage
            %
            %   [mu, sigma] = obj.get_result(i);
            %       gets the result of a block. In the output, mu 
            %       and sigma are mean vector and covariance matrix 
            %       of the i-th component.
            %            
            %   [mu, sigma] = obj.get_result();
            %       gets all results. In the output, mu and sigma are
            %       both K x 1 cell arrays. mu{k} and sigma{k} correspond
            %       to the mean vector and covariance matrix of the 
            %       k-th block.
            %
            
            if nargin < 2
                K = obj.nblocks;
                p_Js = obj.pre_Js;
                p_hs = obj.pre_hs;
                
                mu = cell(K, 1);
                sigma = cell(K, 1);
                for i = 1 : K
                    p_J = p_Js{i};
                    p_h = p_hs{i};
                    
                    C = inv(p_J);
                    mu{i} = C * p_h; %#ok<MINV>
                    sigma{i} = C;
                end
            else
                p_J = obj.pre_Js{i};
                p_h = obj.pre_hs{i};
                
                if nargout <= 1
                    mu = p_J \ p_h;
                else
                    sigma = inv(p_J);
                    mu = sigma * p_h; %#ok<MINV>
                end                
            end
        end
        
        
        
    end
    
    
    
end