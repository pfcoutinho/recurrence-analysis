function [RP, DM] = recplot(m, t, threshold, dtype, varargin)
%RECPLOT Recurrence plot
%   Generates the recurrence plot of a time series or the cross-recurrence plot
%       between two time series.
%
% SYNTAX
%   [RP, DM]   = recplot(de, tau, threshold, dfunc, x)
%   [XRP, XDM] = recplot(de, tau, threshold, dfunc, x, y)
%
% INPUT
%   m           - embedding dimension
%   t           - time delay
%   threshold   - threshold (scalar, or 1 x 2 or 2 x 1 vector)
%   dtype       - distance function
%   x, y        - time series (vectors)
%
% OUTPUT
%   RP / XRP    - recurrence or cross-recurrence plot (logical matrix)
%   DM / XDM    - distance matrix
%
% REFERENCES
%   Recurrence Plots. Site: <http://www.agnld.uni-potsdam.de/~marwan/ 
%       matlab-tutorials/html/rp.html>
%   N. Marwan, M. C. Romano, M. Thiel, and J. Kurths, Recurrence Plots for
%       the Analysis of Complex Systems, Physics Reports, Vol. 438, pp.237-329,
%       2007 
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%   Last update: Jan 22, 2020
% ============================================================================ %

    %
    % Checking the number of input and output arguments
    %
    
    % Input
    narginchk(5, 6);
    
    % Output
    nargoutchk(1, 2);
    
    %
    % Getting the parameters of the recurrence plot and the time series
    %
    
    % Recurrence plot parameters
    parameters.embeddingDimension = m;
    parameters.timeDelay          = t;
    parameters.threshold          = threshold;
    parameters.distance              = dtype;
    
    % Time series
    x = varargin{1};
    if nargin == 6
        y = varargin{2};
    end
    
    %
    % Checking the arguments (parameters of the recurrence plot and the time
    % series)
    %
    if(nargin == 5)
        FLAG_CORRIDOR = checkarguments(parameters, x);
    elseif(nargin == 6)
        FLAG_CORRIDOR = checkarguments(parameters, x, y);
    end
    
    %
    % Obtain the distance matrix
    %
    switch nargin
        case 5  % recurrence plot;
            DM = distmatrix(parameters, x);
        case 6  % cross recurrence plot;
            DM = distmatrix(parameters, x, y);
    end
    
    % Apply the threshold
    if ~FLAG_CORRIDOR
        % thresholded version;
        RP = (abs(DM) <= parameters.threshold);
    else
        % corridor thresholded version;
        RPa = (abs(DM) >= parameters.threshold(1));
        RPb = (abs(DM) <= parameters.threshold(2));
        RP = and(RPa, RPb);
    end
        
end % END of recplot()

function [DM] = distmatrix(parameters, varargin)
    %
    % Getting the time series
    %
    x = varargin{1};
    if(nargin == 3)
        y = varargin{2};
    end

    %
    % Distance matrix for one time series
    %
    if nargin == 2
        M_MAX = 1e4;
        
        % number of state vectors
        M = length(x)-parameters.timeDelay*(parameters.embeddingDimension-1);

        if M <= M_MAX
            % state vectors;
            X = zeros(M, parameters.embeddingDimension);
            for ii = 1:parameters.embeddingDimension
                X(1:end, ii) = x((1:M)+parameters.timeDelay*(ii-1));
            end

            % replicate state vectors;
            X1 = repmat(X, M, 1);
            X2 = reshape(repmat(X(:),1,M)', M*M, parameters.embeddingDimension);

            % compute the similarity or distance between the state vectors;
            switch parameters.distance
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
                    DM = sum(abs(X1-X2), 2)/parameters.embeddingDimension;
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
            DM = reshape(DM, M, M);

        else
            % pre-allocating the RP matrix;
            DM = zeros(M, M);

            for ii = 1:M-1
                x1 = x(ii:parameters.timeDelay:(ii-1)+parameters.embeddingDimension*parameters.timeDelay)';          % state vector #1;
                for jj = ii+1:M
                    x2 = x(jj:parameters.timeDelay:(jj-1)+parameters.embeddingDimension*parameters.timeDelay)' ;     % state vector #2;

                    switch parameters.distance
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
                            DM(ii, jj) = sum(abs(x1-x2), 2)/parameters.embeddingDimension;
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
        M = length(x)-parameters.timeDelay*(parameters.embeddingDimension-1);
        N = length(y)-parameters.timeDelay*(parameters.embeddingDimension-1);

        if (M*N) <= MxN_MAX
            % state vectors;
            X = zeros(M,parameters.embeddingDimension);
            Y = zeros(N,parameters.embeddingDimension);
            for ii = 1:parameters.embeddingDimension
                X(1:end, ii) = x((1:M)+parameters.timeDelay*(ii-1));
                Y(1:end, ii) = y((1:N)+parameters.timeDelay*(ii-1));
            end
            % replicate state vectors;
            XX = reshape(repmat(X(1:end), 1, N)', M*N, parameters.embeddingDimension);
            YY = repmat(Y, M, 1);

            % compute the similarity or distance between the state vectors;
            switch parameters.distance
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
                    DM = sum(abs(XX-YY),2)/parameters.embeddingDimension;
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
            DM = reshape(DM, M, N);

        else
            % pre-allocating the distance matrix;
            DM = zeros(M, N);

            for ii = 1:M
                xx = x(ii:parameters.timeDelay:(ii-1)+parameters.embeddingDimension*parameters.timeDelay)';       % state vector #1;
                for jj = 1:N
                    yy = y(jj:parameters.timeDelay:(jj-1)+parameters.embeddingDimension*parameters.timeDelay)' ;  % state vector #2;
                    % unthresholded version;
                    switch parameters.distance
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
                            DM(ii, jj) = sum(abs(xx-yy), 2)/parameters.embeddingDimension;
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
        
end % END of distmatrix()

function FLAG_CORRIDOR = checkarguments(parameters, varargin)

    x = varargin{1};
    if nargin==3
        y = varargin{2};
    end
    
    %
    % Checking arguments
    %
    
    % Embedding dimension
    if isnumeric(parameters.embeddingDimension)
        if ~isscalar(parameters.embeddingDimension)
            error("Embeddind dimension must be scalar");
        end 
    else
        error("Embedding dimension must be numeric");
    end
    
    if mod(parameters.embeddingDimension,1) ~= 0
        error("Embedding dimension must be a natural number");
    end
    
    if parameters.embeddingDimension < 1
        error("Embedding dimension must be greather than or equal to 1"); 
    end
    
    % Time delay
    if(isnumeric(parameters.timeDelay))
        if(~isscalar(parameters.embeddingDimension))
            error("Delay time must be scalar");
        end 
    else
        error("Delay time must be numeric");
    end

    if(mod(parameters.timeDelay, 1) ~= 0)
        error("Delay time must be a natural number");
    end

    
    if(parameters.embeddingDimension == 1)
        parameters.timeDelay = 1;        % for computational purpose;
    end    
           
    if(parameters.embeddingDimension > 1 && parameters.timeDelay < 1)
        err_str = ['If embedding dimension is greater than 1, ', ...
                    'tau value must be greater than or equal to 1'];
        error(err_str);
    end
    
    if parameters.timeDelay < 0
        error('Delay must be greater than or equal to 0');
    elseif parameters.timeDelay == 0
        parameters.timeDelay = 1;        % for computational purpose;
    end
    
    % Threshold
    if ~isnumeric(parameters.threshold)
        error('Threshold must be numeric');
    end
    
    if ~isempty(parameters.threshold)
        if isnumeric(parameters.threshold)
            if isscalar(parameters.threshold)
                FLAG_CORRIDOR = false;
            else
                % corridor threhsold;
                if (size(parameters.threshold, 1) == 1 && ...
                        size(parameters.threshold,2) == 2) || ...
                        (size(parameters.threshold, 1) == 2 && ...
                        size(parameters.threshold,2) == 1)
                    if size(parameters.threshold,1) == 2
                        parameters.threshold = ...
                            parameters.threshold';
                    end                        
                    FLAG_CORRIDOR = true;
                else
                    error('Threshold must be 1 x 2 or 2 x 1 vector.');
                end
            end
        else
            error('Threshold (epsilon) must be numeric');
        end
    end
    
    % similarity or distance function;
    if ~ischar(parameters.distance)
        err_str = ['Similarity or distance function must be ', ...
            'specified by a string'];
        error(err_str);
    end   
    
    DFUNCS = {'L1','L2','Linf','cos','Gower','Lorentzian', ...
        'intersec','inner_prod','Jaccard','Dice'};
    if ~ismember(parameters.distance,DFUNCS)
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
    
    if size_x(1) < parameters.timeDelay*(parameters.embeddingDimension-1)
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

        if size_y(1) < parameters.timeDelay*(parameters.embeddingDimension-1)
            error('Insufficient number of samples (2nd data series)');
        end
    end
    
end % END of checkarguments()
