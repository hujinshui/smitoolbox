classdef ltnet
    % The class to represent a linear transformation network
    %
    % A linear transformation network is formulated as a multi input
    % multi output linear map, as defined by follows.
    %
    % It comprises m input vectors in a vector space of dimension d_in,
    % and n output vectors in a vector space of dimension d_out.
    % Given the input vectors, and the associate coefficients c_{ij},
    % the j-th output vector y_j is computed by
    % 
    %   y_j = \sum_{i=1}^m c_{ij} A_{ij} x_i.
    %
    % Here i = 1,...,m, and j = 1,...,n. Each matrix A_{ij} is a 
    % d_out x d_in matrix. 
    %
    
    % Created by Dahua Lin, on Apr 4, 2010.
    %
    
    properties(GetAccess='public', SetAccess='private')
        As;         % the array storing all transforms
        d_in;       % the dimension of input space
        d_out;      % the dimension of output space
        n_in;       % the number of input vectors (if [], it is arbitrary)
        n_out;      % the number of output vectors {if [], it is arbitrary)
        
        tlevel;     % the level at that the transform is tied
                    % 0: A_{ij} is different for different i and j
                    % 1: for each i, A_{ij} is the same for all j
                    % 2: for each j, A_{ij} is the same for all i
                    % 3: A_{ij} is all the same

        form;       % the form of As
                    % 'nrm':    the normal form
                    % 'tran':   the transpose form
                    % 'row':    the vector form when A_{ij} is row
                    % 'col':    the vector form when A_{ij} is column
                    % 'sca':    the scalar form     
                    
        preS;       % pre-computed structure
    end
    
    methods
        
        %% net construction
        
        function obj = ltnet(As, form, tlevel, d)
            % Constructs a linear transform network object
            %
            %   obj = ltnet(As, form, tlevel);
            %   obj = ltnet(As, form, tlevel, d);
            %       constructs a linear transform network object.
            %       The size of As depends on the specific form and
            %       the level of tying (tlevel).
            %
            %       The form is a string that can take either of the
            %       following values:
            %       - 'nrm': the normal form. Each A_{ik} is given by
            %         its matrix representation. 
            %       - 'tran': the transpose form. Each A_{ik} is given by
            %         its transposed matrix representation. 
            %       - 'row': the compact vector form. This applies only
            %         when d_out == 1, thus each A_{ik} is essentially
            %         a row vector.
            %       - 'col': the compact vector form. This applies only
            %         when d_in == 1, thus each A_{ik} is essentially
            %         a column vector.
            %       - 'sca': the scalar form. Each A_{ik} is a scalar.
            %
            %       The value of tlevel can be either of 0, 1, 2, and 3.
            %       - 0:  A_{ij} are respectively specified for different
            %             i and j.
            %       - 1:  A_{ij} are the same for each i.
            %       - 2:  A_{ij} are the same for each j.
            %       - 3:  A_{ij} are al the same.
            %
            %       Depending on the form and tlevel, the size of As
            %       is specified as follows:
            %       For form 'nrm':
            %       - tlevel = 0:  d x r x n x m
            %       - tlevel = 1:  d x r x m
            %       - tlevel = 2:  d x r x n
            %       - tlevel = 3:  d x r
            %       For form 'tran':
            %       - tlevel = 0:  r x d x n x m
            %       - tlevel = 1:  r x d x m
            %       - tlevel = 2:  r x d x n
            %       - tlevel = 3:  r x d
            %       For form 'row':
            %       - tlevel = 0:  r x n x m
            %       - tlevel = 1:  r x m
            %       - tlevel = 2:  r x n
            %       - tlevel = 3:  r x 1
            %       For form 'col':
            %       - tlevel = 0:  d x n x m
            %       - tlevel = 1:  d x m
            %       - tlevel = 2:  d x n
            %       - tlevel = 3:  d x 1
            %       For form 'sca':
            %       - tlevel = 0:  n x m
            %       - tlevel = 1:  m x 1
            %       - tlevel = 2:  n x 1
            %       - tlevel = 3:  1 x 1
            %
            %       When the form is a scalar, an additional argument d
            %       is needed to specify both the input and output
            %       dimension (which are equal).
            %
            
            
            %% verify input arguments
            
            assert(isfloat(As) && ndims(As) <= 4, ...
                'ltnet:invalidarg', ...
                'As should be a numeric array with ndims(As) <= 4.');
            
            assert(ischar(form) && size(form, 1) == 1 && ndims(form) == 2, ...
                'ltnet:invalidarg', ...
                'form should be a string.');
            
            assert(isnumeric(tlevel) && isscalar(tlevel) && tlevel == fix(tlevel) ...
                && tlevel >= 0 && tlevel <= 3, ...
                'ltnet:invalidarg', ...
                'tlevel should be an integer scalar in {0,1,2,3}.');
            
            switch form
                case {'nrm', 'tran'}
                    switch tlevel
                        case 0
                            obj.n_out = size(As, 3);
                            obj.n_in = size(As, 4);
                            max_ndim = 4;
                        case 1
                            obj.n_in = size(As, 3);
                            max_ndim = 3;
                        case 2
                            obj.n_out = size(As, 3);
                            max_ndim = 3;
                        case 3
                            max_ndim = 2;
                    end

                    if strcmp(form, 'nrm')
                        obj.d_in = size(As, 2);
                        obj.d_out = size(As, 1);
                    else
                        obj.d_in = size(As, 1);
                        obj.d_out = size(As, 2);
                    end
                    
                case {'row', 'col'}
                    switch tlevel
                        case 0
                            obj.n_out = size(As, 2);
                            obj.n_in = size(As, 3);
                            max_ndim = 3;
                        case 1
                            obj.n_in = size(As, 2);
                            max_ndim = 2;
                        case 2
                            obj.n_out = size(As, 2);
                            max_ndim = 2;
                        case 3
                            max_ndim = 1;
                    end
                    
                    if strcmp(form, 'row')
                        obj.d_in = size(As, 1);
                        obj.d_out = 1;
                    else
                        obj.d_in = 1;
                        obj.d_out = size(As, 1);
                    end
                    
                case 'sca'
                    assert(nargin == 4, 'ltnet:invalidarg', ...
                        'd should be specified when the form is ''sca''.');
                    
                    assert(isnumeric(d) && isscalar(d) && d == fix(d) && d >= 1, ...
                        'ltnet:invalidarg', ...
                        'd should be a positive integer scalar.');
                    
                    switch tlevel
                        case 0
                            obj.n_out = size(As, 1);
                            obj.n_in = size(As, 2);
                            max_ndim = 2;
                        case 1
                            obj.n_in = size(As, 1);
                            max_ndim = 1;
                        case 2
                            obj.n_out = size(As, 1);
                            max_ndim = 1;
                        case 3
                            max_ndim = 0;
                    end
                    
                    obj.d_in = d;
                    obj.d_out = d;
                    
                otherwise
                    error('ltnet:invalidarg', ...
                        'The form %s is unsupported.', form);
            end
            
            if max_ndim >= 2
                assert(ndims(As) <= max_ndim, 'ltnet:invalidarg', ...
                    'As should have ndims(As) <= %d under specified form and tlevel.', max_ndim);
            elseif max_ndim == 1
                assert(ndims(As) == 2 && size(As, 2) == 1, ...
                    'ltnet:invalidarg', ...
                    'As should be a column vector under specified form and tlevel.');
            elseif max_ndim == 0
                assert(isscalar(As), 'ltnet:invalidarg', ...
                    'As should be a scalar under specified form and tlevel.');
            end
            
            obj.As = As;
            obj.form = form;
            obj.tlevel = tlevel;
                              
            % Do necessary pre-computation
            
            tl = obj.tlevel;
            r = obj.d_in;
            d = obj.d_out;
            m = obj.n_in;
            n = obj.n_out;
            
            if tl == 0
                T = [];
                switch form
                    case 'nrm'                        
                        T = reshape(As, [d*r, n*m]);
                    case 'tran'
                        T = permute(As, [2 1 3 4]);
                        T = reshape(T, [d*r, n*m]);
                    case 'row'                       
                        T = reshape(As, [r, n*m]);
                    case 'col'
                        T = reshape(As, [d, n*m]);
                end
                obj.preS = T;                
                                
            elseif tl == 1 || tl == 2              
                if strcmp(form, 'tran')
                    obj.preS = permute(As, [2 1 3]);
                end                
            end
            
        end
        
        
        %% net evaluation
        
        function Y = evaluate(obj, X, C)
            % Evaluates the output vectors from input
            %
            %   Y = obj.evaluate(X, C);
            %       evaluates the output vectors based on the input
            %       vectors given by X and the associate coefficients
            %       given by C.
            %
            %       X should be a numeric matrix of size d_in x m, with
            %       each column being an input vector. 
            %       
            %       C should be a numeric matrix of size m x n, where 
            %       C(i, j) associates the i-th input with the j-th output.
            %
            %       In the output, Y is a numeric matrix of size d_out x n,
            %       where Y(:, j) gives the j-th output vector.
            %       
            %       Note that, when tie level is 1, n can be arbitrary,
            %       when tie level is 0, both m and n can be arbitrary.
            %       However, the constraint that should always be met is
            %       that size(C, 1) == size(X, 2) and 
            %       size(C, 2) == size(Y, 2). 
            %
            
            %% verify input arguments
            
            r = obj.d_in;
            d = obj.d_out;
            m = obj.n_in;
            n = obj.n_out;
            tl = obj.tlevel;
            
            switch tl
                case 0
                    assert(isfloat(C) && isequal(size(C), [m n]), ...
                        'ltnet:evaluate:invalidarg', ...
                        'C should be a numeric matrix of size m x n.');
                    
                    assert(isfloat(X) && isequal(size(X), [r m]), ...
                        'ltnet:evaluate:invalidarg', ...
                        'X should be a numeric matrix of size r x m.');
                case 1
                    assert(isfloat(C) && size(C, 1) == m, ...
                        'ltnet:evaluate:invalidarg', ...
                        'C should be a numeric matrix with m rows.');
                    n = size(C, 2);
                    
                    assert(isfloat(X) && isequal(size(X), [r m]), ...
                        'ltnet:evaluate:invalidarg', ...
                        'X should be a numeric matrix of size r x m.');
                case 2
                    assert(isfloat(C) && size(C, 2) == n, ...
                        'ltnet:evaluate:invalidarg', ...
                        'C should be a numeric matrix with n rows.');
                    m = size(C, 1);
                    
                    assert(isfloat(X) && isequal(size(X), [r m]), ...
                        'ltnet:evaluate:invalidarg', ...
                        'X should be a numeric matrix of size r x m.');
                case 3
                    assert(isfloat(C), ...
                        'ltnet:evaluate:invalidarg', ...
                        'C should be a numeric matrix.');
                    [m, n] = size(C);
                    
                    assert(isfloat(X) && isequal(size(X), [r m]), ...
                        'ltnet:evaluate:invalidarg', ...
                        'X should be a numeric matrix of size r x m.');
            end
                    
            %% do evaluation
            
            A = obj.As;
            frm = obj.form;
            
            if ~strcmp(frm, 'sca')
                if tl == 0
                    T = obj.preS;
                    cs = reshape(C.', [1, m * n]);
                    T = bsxfun(@times, T, cs);
                    T = reshape(T, [d, r, n, m]);
                    T = tilemat(T);
                    
                    Y = T * X(:);
                    Y = reshape(Y, [d, n]);
                    
                elseif tl == 1                      
                    switch frm
                        case 'nrm'                            
                            Z = mmvmult(A, X);
                        case 'tran'
                            Z = mmvmult(obj.preS, X);
                        case 'row'
                            Z = sum(A .* X, 1);
                        case 'col'
                            Z = bsxfun(@times, A, X);
                    end  
                    
                    Y = Z * C;
                    
                elseif tl == 2
                    CX = X * C;
                    
                    switch frm
                        case 'nrm'
                            Y = mmvmult(A, CX);
                        case 'tran'
                            Y = mmvmult(obj.preS, CX);
                        case 'row'
                            Y = sum(A .* CX, 1);
                        case 'col'
                            Y = bsxfun(@times, A, CX);
                    end
                                        
                else % tl == 3
                    switch frm 
                        case {'nrm', 'col'}
                            Z = A * X;
                        case {'tran', 'row'}
                            Z = A' * X;
                    end
                            
                    Y = Z * C;
                end
                
            else  % scalar form
                
                switch tl
                    case 0
                        W = C .* A.';
                    case 1
                        W = bsxfun(@times, C, A);
                    case 2
                        W = bsxfun(@times, C, A.');
                    case 3
                        W = C * A;
                end
                               
                Y = X * W;                
            end
            
        end
        
    end
        
end

