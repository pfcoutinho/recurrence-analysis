classdef Recurrence
%RECURRENCE Recurrence
%   This class contains the recurrence parameters that are used to generate the
%   distance matrix and the recurrence plot, which are the embedding dimension, 
%   the time delay, the threshold and the norm.
%
% SYNTAX
%   obj = Recurrence(embeddingDimension, timeDelay, normType);
%   obj = Recurrence(embeddingDimension, timeDelay, threshold, normType);
%
% PROPERTIES
%   embeddingDimension
%       Embedding dimension (positive integer greater than or equal to 1)
%
%   timeDelay
%       Time delay (positive integer greater than or equal to 1)
%
%   threshold
%       Threshold (non-negative real number or interval)
%
%   normType
%       Norm used to compute the distance between state-space vectors
%
% METHODS (PUBLIC)
%   Recurrence()
%       Class constructor function
% 
%   dm()
%       Self-distance matrix or cross distance matrix
%
%   rp()
%       Recurrence plot or cross recurrence plot
%
%   plotr()
%       Plot distance matrix or recurrence plot
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Feb 15, 2020
% ============================================================================ %

    properties
        % Recurrence parameters
        embeddingDimension  % Embedding dimension
        timeDelay           % Time delay
        threshold           % Threshold
        normType            % Norm (L1, L2, L-infinity)
    end % END properties
    
    methods
        %
        % Constructor method
        %
        
        function obj = Recurrence(varargin)
        %RECURRENCE
        %   To obtain the distance matrix, 3 parameters must be given: embedding
        %   dimension, time delay, and a norm; the recurrence plot needs 4:
        %   embedding dimension, time delay, threshold and a norm.
        % -------------------------------------------------------------------- %
            switch nargin
                case 3  % Distance matrix
                    obj.embeddingDimension = varargin{1};
                    obj.timeDelay          = varargin{2};
                    obj.normType           = varargin{3};
                    
                    % Check parameters
                    obj = argchk(obj, 'distance matrix');
                case 4  % Recurrence plot
                    obj.embeddingDimension = varargin{1};
                    obj.timeDelay          = varargin{2};
                    obj.threshold          = varargin{3};
                    obj.normType           = varargin{4};
                    
                    % Check parameters
                    obj = argchk(obj, 'recurrence plot');
                otherwise
                    ERR_MSG = "Wrong number of input arguments.";
                    error(ERR_MSG);
            end
            
        end % END Recurrence()
        
        %
        % Distance matrix
        %
        
        function D = dm(obj, varargin)
        %DM Distance matrix
        %
        % -------------------------------------------------------------------- %
            % Check recurrence object
            if(~isa(obj, 'Recurrence'))
                ERR_MSG = "First argument must be an object of the class Recurrence.";
                error(ERR_MSG);
            end

            % Validate data
            switch nargin
                case 2
                    obj.chkdata(varargin{1});

                    % Number of state-space vectors
                    N = numel(varargin{1}) - obj.timeDelay ...
                        * (obj.embeddingDimension - 1);
                    if(N < 1)
                        ERR_MSG = "Time series is too short.";
                        error(ERR_MSG);
                    end
                case 3
                    obj.chkdata(varargin{1});
                    obj.chkdata(varargin{2});

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
            
            % Obtain the self-distance / cross distance matrix
            switch nargin
                case 2  % If there's only one time series
                    D = sdm(obj, varargin{1:end});
                case 3  % If there are two time series
                    D = cdm(obj, varargin{1:end});
            end
        end % END dm()
        
        %
        % Recurrence plot
        %

        function RP = rp(obj, varargin)
        %RP Recurrence plot
        %
        % -------------------------------------------------------------------- %
            % Obtain the distance matrix (DM)
            DM = dm(obj, varargin{1:end});
            
            if(obj.threshold(1) == 0)
                RP = (DM <= obj.threshold(2));
            else
                RP = and(DM >= obj.threshold(1), DM <= obj.threshold(2));
            end
        end % END rp()
        
        %
        % Recurrence Quantification Analysis
        %
        
        
        %
        % Sets
        % 
        
        function obj = set.embeddingDimension(obj, newEmbeddingDimension)
        %SET.EMBEDDINGDIMENSION Set the value of the embedding dimension
        %   Changes the value of the embedding dimension.
        % -------------------------------------------------------------------- %
            obj.embeddingDimension = newEmbeddingDimension;
        end % END set.embeddingDimension()
        
        function obj = set.timeDelay(obj, newTimeDelay)
        %SET.TIMEDELAY Set the value of the time delay
        %   Changes the value of the time delay.
        % -------------------------------------------------------------------- %
            obj.timeDelay = newTimeDelay;
        end % END set.timeDelay()
        
        function obj = set.threshold(obj, newThreshold)
        %SET.THRESHOLD Set the value of the threshold
        %   Changes the value of the threshold.
        % -------------------------------------------------------------------- %
            obj.threshold = newThreshold;
        end % END set.threshold()
        
        function obj = set.normType(obj, newNormType)
        %SET.NORMTYPE Set the norm
        %   Changes the norm.
        % -------------------------------------------------------------------- %
            obj.normType = newNormType;
        end % END set.normType()
        
    end % END public methods
    
    methods (Access = private)
        %
        % Checking parameters
        %
        
        function obj = argchk(obj, opt)
        %ARGCHK Check recurrence parameters
        %   This function validates the embedding dimension, the time delay, the
        %   threshold, and the norm.
        % -------------------------------------------------------------------- %
            switch opt
                case 'distance matrix'
                    obj = chkembeddingdimension(obj);
                    obj = chktimedelay(obj);
                    obj = chknorm(obj);
                case 'recurrence plot'
                    obj = chkembeddingdimension(obj);
                    obj = chktimedelay(obj);
                    obj = chkthreshold(obj);
                    obj = chknorm(obj);
            end
        end % END checkparameters()
        
        function obj = chkembeddingdimension(obj)
        %CHKEMBEDDINGDIMENSION
        %   This function validates the embedding dimension. The value must be a
        %   positive integer number greater than or equal to 1.
        % -------------------------------------------------------------------- %
            ERR_MSG = "Invalid parameter: embedding dimension.";
        
            if(isempty(obj.embeddingDimension))
                ERR_MSG = strcat(ERR_MSG, " Argument cannot be empty.");
            end
            
            if(ischar(obj.embeddingDimension))
                obj.embeddingDimension = sprintf("%c", obj.embeddingDimension);
            end
        
            if(isstring(obj.embeddingDimension))
                obj.embeddingDimension = str2double(obj.embeddingDimension);
            end
                
            [m,n] = size(obj.embeddingDimension);
            if(m ~= 1 || n ~= 1)
                ERR_MSG = strcat(ERR_MSG, " It must be a positive integer", ...
                    " number greater than or equal to 1.");
                error(ERR_MSG);
            end
            
            % Embedding dimension < 1?
            if(obj.embeddingDimension < 1)
                ERR_MSG = strcat(ERR_MSG, " It must be a positive integer", ...
                    " number greater than or equal to 1.");
                error(ERR_MSG);
            end
            
            % Integer?
            if(mod(obj.embeddingDimension, 1) ~= 0)
                ERR_MSG = strcat(ERR_MSG, " It must be a positive integer", ...
                    " number.");
                error(ERR_MSG);
            end
        end % END chkembeddingdimension()
        
        function obj = chktimedelay(obj)
        %CHKTIMEDELAY Check time delay
        %   This function validates the time delay, which is integer number
        %   greater than or equal to 0. If the embedding dimension is equal to
        %   1, then the value of the time delay doesn't matter (in this case, we
        %   set time delay to 0).
        %
        %   User can input the time delay value as a number (e.g., 2), as string
        %   (e.g., "2"), or as a char (e.g., '2').
        % -------------------------------------------------------------------- %
            ERR_MSG = "Invalid parameter: time delay.";
            
            if(isempty(obj.timeDelay))
                ERR_MSG = strcat(ERR_MSG, " Argument cannot be empty.");
                error(ERR_MSG);
            end
            
            if(ischar(obj.timeDelay))
                obj.timeDelay = sprintf("%c", obj.timeDelay);
            end
        
            if(isstring(obj.timeDelay))
                obj.timeDelay = str2double(obj.timeDelay);
                if(isempty(obj.timeDelay))
                    error(ERR_MSG);
                end
            end
            
            if(isnumeric(obj.timeDelay))
                % Negative number?
                if(obj.timeDelay < 0)
                    ERR_MSG = strcat(ERR_MSG, " Value must be positive.");
                    error(ERR_MSG);
                end
                
                % Vector?
                [m, n] = size(obj.timeDelay);
                if(m > 1 || n > 1)
                    ERR_MSG = strcat(ERR_MSG, " Value must be a scalar.");
                    error(ERR_MSG);
                end
                
                % Not an integer?
                if(mod(obj.timeDelay, 1) ~= 0)
                    ERR_MSG = strcat(ERR_MSG, " Value must be an integer.");
                    error(ERR_MSG);
                end
            end
            
            if(obj.embeddingDimension == 1 && obj.timeDelay ~= 0)
                obj.timeDelay = 0;
            end
            
            if(obj.embeddingDimension > 1 && obj.timeDelay == 0)
                ERR_MSG = strcat(ERR_MSG, " Value must be greater than", ...
                    " or equal to 1.");
                error(ERR_MSG);
            end
        end % END chktimedelay()
        
        function obj = chkthreshold(obj)
        %CHKTHRESHOLD Check the threshold parameter
        %   Threshold value
        % -------------------------------------------------------------------- %
            ERR_MSG = "Invalid parameter: threshold.";
            
            if(isempty(obj.threshold))
                ERR_MSG = strcat(ERR_MSG, " Argument cannot be empty.");
                error(ERR_MSG);
            end
        
            if(ischar(obj.threshold))
                obj.threshold = sprintf("%c", obj.threshold);
            end
            
            if(isstring(obj.threshold))
                obj.threshold = str2double(obj.threshold);
            end
            
            [m, n] = size(obj.threshold);
            if(m == 1 && n == 1)
                obj.threshold(2) = 0;
                [m, n] = size(obj.threshold);
            end
            
            if(m < 0 || m > 2 || n < 0 || n > 2)
                ERR_MSG = strcat(ERR_MSG, " It must be a single real", ...
                    " number or an interval (e.g., [0.1, 0.25]).");
                error(ERR_MSG);
            end
            
            if(m == 2 && n == 2)
                ERR_MSG = strcat(ERR_MSG, " It must be a single real", ...
                    " number or an interval (e.g., [0.1, 0.25]).");
                error(ERR_MSG);
            end
            
            if((m == 2 && n == 1) || (m == 1 && n == 2))
                % Corridor thresholded version of the recurrence plot
                obj.threshold = sort(obj.threshold, 'ascend');
            end            
        end % END chkthreshold()
        
        function obj = chknorm(obj)
        %CHKNORM Check the norm used to calculate the distance between
        %state-space vectors
        %   This function validates the norm. It must be choosen among L1
        %   (a.k.a. taxicab norm), L2 (a.k.a. Euclidean norm), and L-infinity
        %   (a.k.a. maximum norm).
        %   
        %   User can use these norms as follows:
        %       * as strings: "L1", "L2", and "L3";
        %       * as chars: 'L1', 'L2', and 'L3';
        %       * as numbers: 1, 2, and 3 (respectively to L1, L2, and L3).
        % -------------------------------------------------------------------- %
            ERR_MSG = "Invalid parameter: norm.";
            
            if(isstring(obj.normType))
                obj.normType = char(obj.normType);
            end
            
            if(isnumeric(obj.normType))
                if(mod(obj.normType, 1) == 0)
                    if(obj.normType < 1 || obj.normType > 3)
                        ERR_MSG = strcat(ERR_MSG, " Norm can be an integer ", ...
                            "number between 1 and 3: 1 = L1, 2 = L2, and ", ...
                            "3 = L-infinity.");
                        error(ERR_MSG);
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
                    ERR_MSG = strcat(ERR_MSG, " Norm can be an integer ", ...
                    "number: 1 = L1,2 = L2, 3 = L-infinity.");
                    error(ERR_MSG);
                end
                
            end
            
            VALID_NORMS = {'L1','l1', 'taxicab', ...
                'L2', 'l2', 'euclidean', 'Euclidean', ...
                'L-infinity', 'l-infinity', 'Linf', 'linf', 'Chebychev'};
            if(~ismember(obj.normType, VALID_NORMS))
                ERR_MSG = strcat(ERR_MSG, " Choose among L1, L2, and ", ...
                    "L-infinity.");
                error(ERR_MSG);
            end

            % Standardize
            switch obj.normType
                case {'l1', 'taxicab'}
                    obj.normType = 'L1';
                case {'l2', 'euclidean', 'Euclidean'}
                    obj.normType = 'L2';
                case {'l-infinity', 'Linf', 'linf', 'Chebychev'}
                    obj.normType = 'L-infinity';
            end
        end % END chknorm()
        
        %
        % Calculation of self and cross-distances
        %
        
        function D = sdm(obj, varargin)
        %SDM Self-distance matrix
        % -------------------------------------------------------------------- %
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
                case {'L-infinity', 'l-infinity', 'L-inf', 'l-inf', ...
                        'Maximum', 'maximum'}
                    D = max(abs(U-V), [], 2);
            end

            D = reshape(D, N, N);
        end % END sdm()

        function D = cdm(obj, varargin)
        %CDM Cross-distance matrix
        % -------------------------------------------------------------------- %
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
                case {'L-infinity', 'l-infinity', 'L-inf', 'l-inf', ...
                        'Maximum', 'maximum'}
                    D = max(abs(U-V), [], 2);
            end

            D = reshape(D, N1, N2);
        end % END cdm()
        
    end % END private methods
    
    methods (Static = true)
        %
        % Plot
        %
        
        function plotr(M)
        %PLOTR Plot recurrence object
        %   Plots the self-distance / cross-distance matrix or the recurrence /
        %   cross-recurrence plot. If the object is a recurrence /
        %   cross-recurrence plot, the image is set to black and white.
        % -------------------------------------------------------------------- %
            imagesc(M)
            
            if(islogical(M))
                % Set colors to B&W
                colormap([1, 1, 1; rand(1, 3)])
            else
                colormap('parula')
                colorbar
            end

            [m, n] = size(M);
            if(m == n)  % symmetric matrix
                axis square
            end

            % Axis
            [m, n] = size(M);
            ax = [1, round(n/2), n];
            ay = [1, round(m/2), m];
            set(gca, 'XTick', ax)
            set(gca, 'YTick', ay)
        end % END plotr()
        
    end % END static methods
    
    methods (Access = private, Static = true)
        %
        % Check the data before calculating the distance matrix or the
        % recurrence plot
        %
        
        function chkdata(data)
        %CHKDATA Check data from time series
        %
        % -------------------------------------------------------------------- %
            if(isempty(data))
                ERR_MSG = "Data from the time series cannot be empty.";
                error(ERR_MSG);
            end

            % Is numeric?
            if(~isnumeric(data))
                ERR_MSG = "Time series must be numeric.";
                error(ERR_MSG);
            end
            
            % Is logical?
            if(islogical(data))
                ERR_MSG = "Time series cannot be logical.";
                error(ERR_MSG);
            end

            % Is too large?
            [m, n] = size(data);
            if(m > 10000 || n > 10000)
                ERR_MSG = strcat("Time series has too many samples. ", ...
                    "Consider using a shorter one (with less than 10,000", ...
                    " samples).");
                error(ERR_MSG);
            end
        end % END chkdata()
        
    end % END private, static methods
    
end
