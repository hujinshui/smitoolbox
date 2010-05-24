classdef gaussd_map < smodelest
    % The class for MAP estimation of Gaussian models (with fixed
    % covariance)
    %    
    
    % Created by Dahua Lin, on Apr 13, 2010
    % Modified by Dahua Lin, on Apr 21, 2010
    %   - based on smodelest
    %
        
    properties(GetAccess='public', SetAccess='private')
        dim;            % dimension of vector space
        cov;            % the observation covariance
        icov;           % the inverse covariance of observation
    end
    
    methods
        
        function obj = gaussd_map(gpri, cov)
            % Constructs a Gaussian MAP estimator
            %
            %   obj = gaussd_map(gpri, cov);
            %       constructs a Gaussian MAP estimator with the 
            %       Gaussian prior given by gpri.
            %
            
            assert((isa(gpri, 'gaussd_cp') || isa(gpri, 'gaussd_mp')) && ...
                gpri.nmodels == 1, ...
                'gaussd_map:invalidarg', ...
                'gpri should be a Gaussian model object with single model.');
            
            if isa(gpri, 'gaussd_mp')
                gpri = to_cp(gpri);
            end
            
            assert(isa(cov, 'pdmat') && cov.num == 1, ...
                'gaussd_map:invalidarg', ...
                'cov should be an object of pdmat with single matrix.');
            
            assert(cov.dim == gpri.dim, ...
                'gaussd_map:invalidarg', ...
                'The dimension of gpri and cov is inconsistent.');
            
            obj.prior = gpri;
            obj.dim = gpri.dim;            
            obj.cov = cov;
            obj.icov = inv(cov);        
        end
        
        
        function Pb = accept(obj, X)
            % Checks the input observation and forms an estimation problem
            %                        
                        
            Pb = gaussd_map_pb(obj, X);        
        end
                
    end
    
end