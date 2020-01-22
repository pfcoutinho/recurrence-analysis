function [DM] = distmatrix(varargin)
%DISTMATRIX Distance matrix
%
% SYNTAX:
%   DM = distmatrix(EMBEDDING_DIMENSION, TIME_DELAY, DTYPE, x)
%   DM = distmatrix(EMBEDDING_DIMENSION, TIME_DELAY, DTYPE, x, y);
%
% INPUT:
%   x, y  - time series (vectors)
%   EMBEDDING_DIMENSION - embedding dimension
%   TIME_DELAY          - time delay
%   DTYPE               - distance function (Linf, L1, L2)
%
% OUTPUT:
%   DM - distance matrix
%
% DEPENDENCIES:
%   None.
%
% CONTACT:
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%   Last update: Jan 22, 2020
% ============================================================================ %

    %
    % Checking the number of input and output arguments
    %
    
    % Number of input arguments must be 4 or 5
    narginchk(4, 5);
    
    % Number of output arguments is 1
    nargoutchk(1, 1);
    % ------------------------------------------------------------------------ %
    
    %
    % Getting the parameters of the distance matrix and the time series
    %
    
    % Distance matrix:
    EMBEDDING_DIMENSION = varargin{1};
    TIME_DELAY          = varargin{2};
    DTYPE               = varargin{3};
    
    % Time series:
    x = varargin{4};
    if(nargin == 5)
        y = varargin{5};
    end
    % ------------------------------------------------------------------------ %
    
    %
    % Checking arguments
    %
    
    % Embedding dimension
    if isnumeric(EMBEDDING_DIMENSION)
        if ~isscalar(EMBEDDING_DIMENSION)
            error("Embeddind dimension must be scalar");
        end 
    else
        error("Embedding dimension must be numeric");
    end
    
    if mod(EMBEDDING_DIMENSION,1) ~= 0
        error("Embedding dimension must be a natural number");
    end
    
    if EMBEDDING_DIMENSION < 1
        error("Embedding dimension must be greather than or equal to 1"); 
    end
    
    % Time delay
    if(isnumeric(TIME_DELAY))
        if(~isscalar(EMBEDDING_DIMENSION))
            error("Delay time must be scalar");
        end 
    else
        error("Delay time must be numeric");
    end

    if(mod(TIME_DELAY, 1) ~= 0)
        error("Delay time must be a natural number");
    end

    
    if(EMBEDDING_DIMENSION == 1)
        TIME_DELAY = 1;        % for computational purpose;
    end    
           
    if(EMBEDDING_DIMENSION > 1 && TIME_DELAY < 1)
        err_str = ['If embedding dimension is greater than 1, ', ...
                    'tau value must be greater than or equal to 1'];
        error(err_str);
    end
    
    if TIME_DELAY < 0
        error('Delay must be greater than or equal to 0');
    elseif TIME_DELAY == 0
        TIME_DELAY = 1;        % for computational purpose;
    end
    
    % similarity or distance function;
    if ~ischar(DTYPE)
        err_str = ['Similarity or distance function must be ', ...
            'specified by a string'];
        error(err_str);
    end   
    
    DFUNCS = {'L1','L2','Linf','cos','Gower','Lorentzian', ...
        'intersec','inner_prod','Jaccard','Dice'};
    if ~ismember(DTYPE,DFUNCS)
        error('Invalid similiarity or distance function');
    end
    
    % time series;
    if ~isnumeric(x)
        error('Data series must be numeric vectors');
    end
    
    size_x = size(x);
    if (size_x(1) ~= 1 && size_x(2) ~= 1)
        error('Data series must be N x 1 (or 1 x N) vector');
    end
    
    if size_x(1) < size_x(2)
        x = x'; % row vector -> column vector;
        size_x = size(x);
    end
    
    if size_x(1) < TIME_DELAY*(EMBEDDING_DIMENSION-1)
        error('Insufficient number of samples (1st data series)');
    end
    
    if nargin == 5
        if ~isnumeric(y)
            error('Data series must be numeric vectors');
        end
        
        size_y = size(y);
        if (size_y(1) ~= 1 && size_y(2) ~= 1)
            error('Data series must be N x 1 (or 1 x N) vector');
        end

        if size_y(1) < size_y(2)
            y = y'; % row vector -> column vector;
            size_y = size(y);
        end

        if size_y(1) < TIME_DELAY*(EMBEDDING_DIMENSION-1)
            error('Insufficient number of samples (2nd data series)');
        end
    end
    
% ============================================================================ %
    %
    % Distance matrix for one time series
    %
    if nargin == 4
        M_MAX = 1e4;
            
        % number of state vectors
        M = length(x)-TIME_DELAY*(EMBEDDING_DIMENSION-1);

        if M <= M_MAX
            % state vectors;
            X = zeros(M, EMBEDDING_DIMENSION);
            for ii = 1:EMBEDDING_DIMENSION
                X(:, ii) = x((1:M)+TIME_DELAY*(ii-1));
            end

            % replicate state vectors;
            X1 = repmat(X, M, 1);
            X2 = reshape(repmat(X(:),1,M)', M*M, EMBEDDING_DIMENSION);

            % compute the similarity or distance between the state vectors;
            switch DTYPE
                case 'Linf'         % Chebyshev L-infinity;
                    DM = max(abs(X1-X2), [], 2);
                case 'L1'           % City block L1;
                    DM = sum(abs(X1-X2), 2);
%                     DM = sum((X1-X2), 2);
                case 'L2'           % Euclidean L2;
                    DM = sqrt(sum((X1-X2).^2, 2));
                case 'cos'          % cosine similarity;
                    DM = sum(X1.*X2,2)./(sqrt(sum(X1.^2, 2)).* ...
                        sqrt(sum(X2.^2, 2)));
                case 'Gower'        % Gower;
                    DM = sum(abs(X1-X2), 2)/EMBEDDING_DIMENSION;
                case 'Lorentzian'   % Lorentzian;
                    DM = sum(log10(1+abs(X1-X2)), 2);
                case 'intersec'     % intersection;
                    DM = sum(abs(X1-X2), 2)/2;
                case 'inner_prod'   % inner product;
                    DM = sum(X1.*X2, 2);
                case 'Jaccard'      % Jaccard distance;
                    DM = sum((X1-X2).^2, 2)./ ...
                        (sum(X1.^2,2)+sum(X2.^2, 2)-sum(X1.*X2, 2));
                case 'Dice'         % Dice distance;
                    DM = sum((X1-X2).^2, 2)./(sum(X1.^2,2)+sum(X2.^2, 2));
            end

            % reshape;
            DM = reshape(DM,M,M);

        else
            % pre-allocating the RP matrix;
            DM = zeros(M, M);

            for ii = 1:M-1
                x1 = x(ii:TIME_DELAY:(ii-1)+EMBEDDING_DIMENSION*TIME_DELAY)';          % state vector #1;
                for jj = ii+1:M
                    x2 = x(jj:TIME_DELAY:(jj-1)+EMBEDDING_DIMENSION*TIME_DELAY)' ;     % state vector #2;

                    switch DTYPE
                        case 'Linf' % L-infinity-norm;
                            DM(ii, jj) = max(abs(x1-x2));
                        case 'L1'   % L1-norm;
                            DM(ii, jj) = sum(abs(x1-x2));
                        case 'L2'   % L2-norm (Euclidean norm);
                            DM(ii, jj) = sqrt(sum((x1-x2).^2));
                        case 'cos'  % cosine similarity;
                            DM(ii, jj) = sum(x1.*x2)/ ...
                                (sqrt(sum(x1.^2))*sqrt(sum(x2.^2)));
                        case 'Gower'        % Gower;
                            DM(ii, jj) = sum(abs(x1-x2), 2)/EMBEDDING_DIMENSION;
                        case 'Lorentzian'   % Lorentzian;
                            DM(ii, jj) = sum(log10(1+abs(x1-x2)), 2);
                        case 'intersec'     % intersection;
                            DM(ii, jj) = sum(abs(x1-x2), 2)/2;
                        case 'inner_prod'   % inner product;
                            DM(ii, jj) = sum(x1.*x2, 2);
                        case 'Jaccard'      % Jaccard distance;
                            DM(ii, jj) = sum((x1-x2).^2, 2)./ ...
                                (sum(x1.^2,2)+sum(x2.^2, 2)-sum(x1.*x2, 2));
                        case 'Dice'         % Dice distance;
                            DM(ii, jj) = sum((x1-x2).^2, 2)./ ...
                                (sum(x1.^2, 2)+sum(x2.^2, 2));
                    end

                    DM(jj, ii) = DM(ii, jj);
                end
            end
        end
        
    else
        MxN_MAX = 1e8;
        
        % number of state vectors;
        M = length(x)-TIME_DELAY*(EMBEDDING_DIMENSION-1);
        N = length(y)-TIME_DELAY*(EMBEDDING_DIMENSION-1);

        if (M*N) <= MxN_MAX
            % state vectors;
            X = zeros(M,EMBEDDING_DIMENSION);
            Y = zeros(N,EMBEDDING_DIMENSION);
            for ii = 1:EMBEDDING_DIMENSION
                X(1:end, ii) = x((1:M)+TIME_DELAY*(ii-1));
                Y(1:end, ii) = y((1:N)+TIME_DELAY*(ii-1));
            end
            % replicate state vectors;
            XX = reshape(repmat(X(1:end), 1, N)', M*N, EMBEDDING_DIMENSION);
            YY = repmat(Y, M, 1);

            % compute the similarity or distance between the state vectors;
            switch DTYPE
                case 'Linf' % l-infinity-norm;
                    DM = max(abs(XX-YY),[],2);
                case 'L1'   % l1-norm;
                    DM = sum(abs(XX-YY),2);
                case 'L2'   % l2-norm (Euclidean norm);
                    DM = sqrt(sum((XX-YY).^2,2));
                case 'cos'  % cosine similarity;
                    DM = sum(XX.*YY,2)./(sqrt(sum(XX.^2,2)).* ...
                        sqrt(sum(YY.^2,2)));
                                case 'Gower'        % Gower;
                    DM = sum(abs(XX-YY),2)/EMBEDDING_DIMENSION;
                case 'Lorentzian'   % Lorentzian;
                    DM = sum(log10(1+abs(XX-YY)),2);
                case 'intersec'     % intersection;
                    DM = sum(abs(XX-YY),2)/2;
                case 'inner_prod'   % inner product;
                    DM = sum(XX.*YY,2);
                case 'Jaccard'      % Jaccard distance;
                    DM = sum((XX-YY).^2,2)./ ...
                        (sum(XX.^2,2)+sum(YY.^2,2)-sum(XX.*YY,2));
                case 'Dice'         % Dice distance;
                    DM = sum((XX-YY).^2,2)./(sum(XX.^2,2)+sum(YY.^2,2));
            end

            % recurrence plot;
            DM = reshape(DM,M,N);

        else
            % pre-allocating the distance matrix;
            DM = zeros(M,N);

            for ii = 1:M
                xx = x(ii:TIME_DELAY:(ii-1)+EMBEDDING_DIMENSION*TIME_DELAY)';       % state vector #1;
                for jj = 1:N
                    yy = y(jj:TIME_DELAY:(jj-1)+EMBEDDING_DIMENSION*TIME_DELAY)' ;  % state vector #2;
                    % unthresholded version;
                    switch DTYPE
                        case 'Linf' % l-infinity-norm;
                            DM(ii, jj) = max(abs(xx-yy));
                        case 'L1'   % l1-norm;
                            DM(ii, jj) = sum(abs(xx-yy));
                        case 'L2'   % l2-norm (Euclidean norm);
                            DM(ii, jj) = sqrt(sum((xx-yy).^2));
                        case 'cos'  % cosine similarity;
                            DM(ii, jj) = sum(xx.*yy,2)./(sqrt(sum(xx.^2,2)).* ...
                                sqrt(sum(xx.^2,2)));
                        case 'Gower'        % Gower;
                            DM(ii, jj) = sum(abs(xx-yy), 2)/EMBEDDING_DIMENSION;
                        case 'Lorentzian'   % Lorentzian;
                            DM(ii, jj) = sum(log10(1+abs(xx-yy)), 2);
                        case 'intersec'     % intersection;
                            DM(ii, jj) = sum(abs(xx-yy), 2)/2;
                        case 'inner_prod'   % inner product;
                            DM(ii, jj) = sum(xx.*yy, 2);
                        case 'Jaccard'      % Jaccard distance;
                            DM(ii, jj) = sum((xx-yy).^2,2)./ ...
                                (sum(xx.^2, 2)+sum(yy.^2,2)-sum(xx.*yy, 2));
                        case 'Dice'         % Dice distance;
                            DM(ii, jj) = sum((xx-yy).^2, 2)./ ...
                                (sum(xx.^2,2)+sum(yy.^2, 2));
                    end
                    
                    DM(jj,ii) = DM(ii, jj);
                end
            end
        end
        
    end
        
end



