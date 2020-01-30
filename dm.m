function D = dm(obj, varargin)
%DM Distance matrix
%   Calculate the distances between state-space vectors
%
% SYNTAX
%   
% INPUT
%   Robj         - Recurrence object
%   timeSeries   - time series
%
% OUTPUT
%   DM          - self-distance or the cross-distance matrix
%
% REFERENCES
%   
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 30, 2020
% ============================================================================ %
    switch nargin
        case 2
            MAX_N_SAMPLES = 1e4;

            % Number of state vectors
            N = length(varargin{1}) ...
                    - obj.timeDelay*(obj.embeddingDimension - 1);

            if N <= MAX_N_SAMPLES
                D = ddc(obj, varargin{1});
            else
                D = idc(obj, varargin{1});
            end
        case 3
            MAX_N_SAMPLES = 1e4;

            % Number of state vectors
            N = length(varargin{1}) ...
                    - obj.timeDelay*(obj.embeddingDimension - 1);

            if N <= MAX_N_SAMPLES
                D = ddc(obj, varargin{1:end});
            else
                D = idc(obj, varargin{1:end});
            end
    end
end % END dm()
        
function D = ddc(obj, varargin)
%DDC Direct distance calculation
%   Generates all state-space vectors and calculate the distance from
%   one to another. Uses a lot of memory (batch computation).
% -------------------------------------------------------------------- %
    switch nargin
        case 2  % There's only one time series
            distanceMatrixType = 'self-distance matrix';
        case 3  % There are two time series
            distanceMatrixType = 'cross-distance matrix';
    end

    switch distanceMatrixType
        case 'self-distance matrix'
            timeSeries = varargin{1};

            % Number of state-space vectors
            N = numel(timeSeries) ...
                    - obj.timeDelay*(obj.embeddingDimension - 1);

            % State-space vectors (each row of X is a state-space
            % vector)
            X = zeros(N, obj.embeddingDimension);
            for i = 1:obj.embeddingDimension
                X(1:end, i) = timeSeries((1:N) + obj.timeDelay*(i - 1));
            end

            % Replicate state-space vectors
            U = repmat(X, N, 1);
            V = reshape(repmat(X, 1, ...
                    N)', obj.embeddingDimension, N*N)';

            % Compute the distance between the state-space vectors
            switch obj.normType
                case {'L1', 'l1'}
                    D = sum(abs(U-V), 2);
                case {'L2', 'l2', 'Euclidean', 'euclidean'}
                    D = sqrt(sum((U-V).^2, 2));
                case {'L-infinity', 'l-infinity', 'L-inf', 'l-inf', ...
                        'Maximum', 'maximum'}
                    D = max(abs(U-V), [], 2);
            end

            D = reshape(D, N, N);

        case 'cross-distance matrix'
            timeSeries1 = varargin{1};
            timeSeries2 = varargin{2};

            % Number of state-space vectors for each time series
            N1 = numel(timeSeries1) ...
                    - obj.timeDelay*(obj.embeddingDimension - 1);
            N2 = numel(timeSeries2) ...
                    - obj.timeDelay*(obj.embeddingDimension - 1);

            % State-space vectors (each row of X and Y is a state-space
            % vector)
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
                case 'L1'
                    D = sum(abs(U-V), 2);
                case 'L2'
                    D = sqrt(sum((U-V).^2, 2));
                case 'L-infinity'
                    D = max(abs(U-V), [], 2);
            end

            D = reshape(D, N1, N2);
    end

end % END ddc()

function D = idc(Robj, varargin)
%IDC Iteractive distance calculation
%   Iteractive process to calculate the distance between state-space
%   vectors. It can take a more time than the ddc() function, but needs
%   less memory.
% -------------------------------------------------------------------- %
    switch nargin
        case 2
            x = varargin{1};
            N = length(x)-Robj.timeDelay*(Robj.embeddingDimension-1);
        case 3
            x = varargin{1};
            y = varargin{2};
    end

    % Matrix pre-allocation
    D = zeros(N, N);

    for i = 1:N-1
        % state vector #1
        x1 = x(i:Robj.timeDelay:(i-1) + Robj.embeddingDimension ...
            * Robj.timeDelay)';
        for j = i+1:N
            % state vector #2
            x2 = x(j:Robj.timeDelay:(j-1) + Robj.embeddingDimension ...
                * Robj.timeDelay)';

            switch Robj.normType
                case 'L1'
                    D(i, j) = sum(abs(x1-x2));
                case 'L2'
                    D(i, j) = sqrt(sum((x1-x2).^2));
                case 'L-infinity'
                    D(i, j) = max(abs(x1-x2));
            end

            D(j, i) = D(i, j);
        end
    end
end % END idc()
