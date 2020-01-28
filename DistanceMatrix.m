classdef DistanceMatrix
%DISTANCEMATRIX Distance matrix
%   This class generates the distance matrix based on the samples of a single
%   time series or between the samples of two time-series.
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 28, 2020
% ============================================================================ %

    properties
        embeddingDimension  % Embedding dimension
        timeDelay           % Time delay
        normType            % Norm (L1, L2, L-infinity)
        M                   % Matrix
    end
    
    methods
        %
        % Class constructor
        %
        
        function obj = DistanceMatrix(embeddingDimension, timeDelay, ...
                normType, varargin)
        
            obj.embeddingDimension = embeddingDimension;
            obj.timeDelay          = timeDelay;
            obj.normType           = normType;
            
            obj.M = dm(obj, varargin{1:end});
        end % END DistanceMatrix()
        
        %
        % Parameters validation
        %
        
        function obj = validateparameter(obj, name)
            switch name
                case 'embedding dimension'
                    validateembeddingdimension(obj);
                case 'time delay'
                    validatetimedelay(obj);
                case 'norm'
                    obj = validatenorm(obj);
            end
        end
        
        function validateembeddingdimension(obj)
        %VALIDATEEMBEDDINGDIMENSION
        % -------------------------------------------------------------------- %
        end % END validateembeddingdimension()
        
        function validatetimedelay(obj)
        %VALIDATETIMEDELAY
        % -------------------------------------------------------------------- %
        end % END validatetimedelay()

        function obj = validatenorm(obj)
        %VALIDATENORM
        % -------------------------------------------------------------------- %   
            if(isstring(obj.normType))
                obj.normType = char(obj.normType);
            elseif(isnumeric(obj.normType))
                if(isscalar(obj.normType))
                    if(obj.normType < 1 && obj.normType > 3)
                        msg = strcat("Norm can be an integer number ", ...
                            "between 1 and 3: 'L1' = 1, 'L2' = 2, and ", ...
                            "'L-infinity' = 3.");
                        error(msg);
                    else
                        switch obj.normType
                            case 1
                                obj.normType = 'L1';
                            case 2
                                obj.normType = 'L2';
                            case 3
                                obj.normType = 'L-infinity';
                        end
                    end
                else
                    msg = strcat("Norm can be an integer number: 'L1' = 1,", ...
                        " 'L2' = 2, or 'L-infinity' = 3.");
                    error(msg);
                end
            end
            
            NORM_TYPES = {'L1','l1','L2','l2','L-infinity','l-infinity'};
            if(~ismember(obj.normType, NORM_TYPES))
                msg = strcat("Invalid norm. Choose among 'L1', 'L2', and ", ...
                    "'L-infinity'.");
                error(msg);
            end

            switch obj.normType
                case 'l1'
                    obj.normType = 'L1';
                case 'l2'
                    obj.normType = 'L2';
                case 'l-infinity'
                    obj.normType = 'L-infinity';
            end
        end % END validatenorm()
        
        %
        % Data validation
        %
        
        function validatedata()
            
        end % END validatedata()
        
        %
        % Distance calculation
        %
        
        function D = dm(obj, varargin)
        %DM Distance matrix
        %   This function calculates the distance matrix.
        % -------------------------------------------------------------------- %
            switch nargin
                case 2
                    MAX_N_SAMPLES = 1e4;
                    
                    % Number of state vectors
                    N = length(varargin{1}) - obj.timeDelay ...
                        * (obj.embeddingDimension - 1);
                    
                    if N <= MAX_N_SAMPLES
                        D = ddc(obj, varargin{1});
                    else
                        D = idc(obj, varargin{1});
                    end
                case 3
            end
        end % END distancematrix()
        
        function D = ddc(obj, varargin)
        %DDC Direct distance calculation
        % -------------------------------------------------------------------- %
            switch nargin
                case 2
                    x = varargin{1};
                    N = length(x)-obj.timeDelay*(obj.embeddingDimension-1);
                case 3
                    x = varargin{1};
                    y = varargin{2};
            end
        
            % State vectors
            X = zeros(N, obj.embeddingDimension);
            for i = 1:obj.embeddingDimension
                X(1:end, i) = x((1:N)+obj.timeDelay*(i-1));
            end

            % Replicate state vectors
            X1 = repmat(X, N, 1);
            X2 = reshape(repmat(X(:),1,N)', N*N, obj.embeddingDimension);

            % Compute the distance between the state vectors
            switch obj.normType
                case 'L1'
                    D = sum(abs(X1-X2), 2);
                case 'L2'
                    D = sqrt(sum((X1-X2).^2, 2));
                case 'L-infinity'
                    D = max(abs(X1-X2), [], 2);
            end
            
            D = reshape(D, N, N); 
        end % END ddc()
        
        function D = idc(obj, x, N)
        %IDC Iteractive distance calculation
        % -------------------------------------------------------------------- %
            switch nargin
                case 2
                    x = varargin{1};
                    N = length(x)-obj.timeDelay*(obj.embeddingDimension-1);
                case 3
                    x = varargin{1};
                    y = varargin{2};
            end
        
            % Matrix pre-allocation
            D = zeros(N, N);

            for i = 1:N-1
                % state vector #1
                x1 = x(i:obj.timeDelay:(i-1) + obj.embeddingDimension ...
                    * obj.timeDelay)';
                for j = i+1:N
                    % state vector #2
                    x2 = x(j:obj.timeDelay:(j-1) + obj.embeddingDimension ...
                        * obj.timeDelay)';

                    switch obj.normType
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
        
        %
        % Plot
        %
        
        function plot(obj)
        %PLOT Plot the distance matrix
        % -------------------------------------------------------------------- %
            
            imagesc(obj.M)
        
            if(~islogical(obj.M))
                colormap('parula')
                colorbar
            end

            if issymmetric(double(obj.M))
                axis square
            end

            % Removes labels
            set(gca, 'XTickLabel', [], 'YTickLabel', [])
        end % END plot()
        
        %
        % Gets and Sets
        % 
        
        function obj = setembeddingdimension(obj, newEmbeddingDimension)
            obj.embeddingDimension = newEmbeddingDimension;
        end % END setembeddingdimension()
        
        
    end
end

