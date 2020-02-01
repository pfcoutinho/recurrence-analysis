function D = dm(obj, varargin)
%DM Distance matrix
%   Calculate the distances between state-space vectors
%
% SYNTAX
%   D = dm(obj, timeSeries)
%   D = dm(obj, timeSeries1, timeSeries2)
%   
% INPUT
%   obj          - Recurrence object
%   timeSeries   - time series
%
% OUTPUT
%   DM           - distance matrix
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 31, 2020
% ============================================================================ %

    % Checking arguments
    argchk(obj, varargin{1:end});
    
    switch nargin
        case 2  % If there's only one time series
            D = sdm(obj, varargin{1:end});
        case 3  % If there are two time series
            D = cdm(obj, varargin{1:end});
    end
end % END dm()

function argchk(obj, varargin)
%ARGCHK Check the arguments of the dm() function
% ============================================================================ %
    % Recurrence object
    if(~isa(obj, 'Recurrence'))
        ERR_MSG = "First argument must be an object of the class Recurrence.";
        error(ERR_MSG);
    end
    
    % Data
    switch nargin
        case 2
            chkdata(varargin{1});
            
            % Number of state-space vectors
            N = numel(varargin{1}) - obj.timeDelay ...
                * (obj.embeddingDimension - 1);
            if(N < 1)
                ERR_MSG = "Time series is too short.";
                error(ERR_MSG);
            end
        case 3
            chkdata(varargin{1});
            chkdata(varargin{2});
            
            % Number of state-space vectors for each time series
            N1 = numel(varargin{1}) - obj.timeDelay ...
                * (obj.embeddingDimension - 1);
            N2 = numel(varargin{1}) - obj.timeDelay ...
                * (obj.embeddingDimension - 1);
            if(N1 < 1 || N2 < 1)
                ERR_MSG = "Time series is too short.";
                error(ERR_MSG);
            end
    end
end

function chkdata(data)
%CHKDATA Check data from time series
%
% ============================================================================ %
    if(isempty(data))
        ERR_MSG = "Data from the time series cannot be empty.";
        error(ERR_MSG);
    end

    if(~isnumeric(data))
        ERR_MSG = "Time series must be numeric.";
        error(ERR_MSG);
    end
    
    [m, n] = size(data);
    if(m > 1e4 || n > 1e4)
        ERR_MSG = strcat("Time series has too many samples. Consider", ...
            "using a shorter one (with less than 10,000 samples).");
        error(ERR_MSG);
    end
    
    
end

function D = sdm(obj, varargin)
%SDM Self-distance matrix
% ============================================================================ %
    timeSeries = varargin{1};

    % Number of state-space vectors
    N = numel(timeSeries) - obj.timeDelay*(obj.embeddingDimension - 1);

    % State-space vectors (each row of X is a state-space
    % vector)
    X = zeros(N, obj.embeddingDimension);
    for i = 1:obj.embeddingDimension
        X(1:end, i) = timeSeries((1:N) + obj.timeDelay*(i - 1));
    end

    % Replicate state-space vectors
    U = repmat(X, N, 1);
    V = reshape(repmat(X, 1, N)', obj.embeddingDimension, N*N)';

    % Compute the distance between the state-space vectors
    switch obj.normType
        case {'L1', 'l1', 'taxicab'}
            D = sum(abs(U-V), 2);
        case {'L2', 'l2', 'Euclidean', 'euclidean'}
            D = sqrt(sum((U-V).^2, 2));
        case {'L-infinity', 'l-infinity', 'L-inf', 'l-inf', 'Maximum', ...
                'maximum'}
            D = max(abs(U-V), [], 2);
    end

    D = reshape(D, N, N);
end % END sdm()

function D = cdm(obj, varargin)
%CDM Cross-distance matrix
% ============================================================================ %
    timeSeries1 = varargin{1};
    timeSeries2 = varargin{2};

    % Number of state-space vectors for each time series
    N1 = numel(timeSeries1) - obj.timeDelay*(obj.embeddingDimension - 1);
    N2 = numel(timeSeries2) - obj.timeDelay*(obj.embeddingDimension - 1);

    % State-space vectors (each row of X and Y is a state-space vector)
    X = zeros(N1, obj.embeddingDimension);
    Y = zeros(N2, obj.embeddingDimension);
    for i = 1:obj.embeddingDimension
        X(1:end, i) = timeSeries1((1:N1) + obj.timeDelay*(i - 1));
        Y(1:end, i) = timeSeries2((1:N2) + obj.timeDelay*(i - 1));
    end

    % Replicate state-space vectors
    U = reshape(repmat(X, 1, N2)', obj.embeddingDimension, N1*N2)';
    V = repmat(Y, N1, 1);

    % Compute the distance between the state-space vectors
    switch obj.normType
        case {'L1', 'l1', 'taxicab'}
            D = sum(abs(U-V), 2);
        case {'L2', 'l2', 'Euclidean', 'euclidean'}
            D = sqrt(sum((U-V).^2, 2));
        case {'L-infinity', 'l-infinity', 'L-inf', 'l-inf', 'Maximum', ...
                'maximum'}
            D = max(abs(U-V), [], 2);
    end

    D = reshape(D, N1, N2);
end % END cdm()
