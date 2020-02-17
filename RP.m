classdef RP < RQA & RNA
%RECURRENCE Recurrence class
%
% SYNTAX
%   obj = RP(embeddingDimension, timeDelay, threshold, normType, data);
%   obj = RP(embeddingDimension, timeDelay, threshold, normType, data1, data2);
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
%   M
%       The recurrence matrix, which is the recurrence plot or the
%       cross-recurrence plot
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
    
    properties (Access = private)
        data                % Data from time series
    end
    
    properties (Dependent = true)
        M                   % Recurrence matrix
    end
    
    methods
        %
        % Class constructor function
        %
        
        function obj = RP(varargin)
        %RECURRENCE Recurrence class constructor
        % -------------------------------------------------------------------- %
            if(nargin < 5 || nargin > 6)
                ERR_MSG = "Wrong number of input arguments.";
                error(ERR_MSG);
            end
        
            % Recurrence parameters
            obj.embeddingDimension = varargin{1};
            obj.timeDelay          = varargin{2};
            obj.threshold          = varargin{3};
            obj.normType           = varargin{4};

            % Data
            obj.data = varargin(5:end);
        end % END Recurrence() 
        
        %
        % Plot
        %
        
        function plotr(obj)
        %PLOTR Plot recurrence object
        %   Plots the self-distance / cross-distance matrix or the recurrence /
        %   cross-recurrence plot. If the object is a recurrence /
        %   cross-recurrence plot, the image is set to black and white.
        % -------------------------------------------------------------------- %
            imagesc(obj.M)
            
            if(islogical(obj.M))
                % Set colors to B&W
                colormap([1, 1, 1; rand(1, 3)])
            else
                colormap('parula')
                colorbar
            end

            [m, n] = size(obj.M);
            if(m == n)  % symmetric matrix
                axis square
            end

            % Axis
            [m, n] = size(obj.M);
            ax = [1, round(n/2), n];
            ay = [1, round(m/2), m];
            set(gca, 'XTick', ax)
            set(gca, 'YTick', ay)
        end % END plotr()
        
        %
        % Recurrence matrix
        %
        
        function M = get.M(obj)
            % Obtain the distance matrix (DM)
            DM = dm(obj);
            
            if(obj.threshold(1) == 0)
                M = (DM <= obj.threshold(2));
            else
                M = and(DM >= obj.threshold(1), DM <= obj.threshold(2));
            end
        end % END get.M()
        
        %
        % Sets
        %
        
        function obj = set.embeddingDimension(obj, value)
        %SET.EMBEDDINGDIMENSION Set the value of the embedding dimension
        %   Validation: embedding dimension must be a integer greater than or 
        %   equal to 1.
        % -------------------------------------------------------------------- %
            obj.embeddingDimension = value;
            
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
                ERR_MSG = strcat(ERR_MSG, " It must be an integer", ...
                    " number greater than or equal to 1.");
                error(ERR_MSG);
            end
            
            % Integer?
            if(mod(obj.embeddingDimension, 1) ~= 0)
                ERR_MSG = strcat(ERR_MSG, " It must be an integer number");
                error(ERR_MSG);
            end
        end % END set.embeddingDimension()
        
        function obj = set.timeDelay(obj, value)
        %SET.TIMEDELAY Set the value of the time delay
        %   Validation: the value of the time delay must be greater than or
        %   equal to 1.
        % -------------------------------------------------------------------- %
            obj.timeDelay = value;
            
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
                if(obj.timeDelay < 1)
                    ERR_MSG = strcat(ERR_MSG, " Value must be greater than 1.");
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

            %{
            if(obj.embeddingDimension == 1 && obj.timeDelay ~= 0)
                obj.timeDelay = 0;
            end
            
            if(obj.embeddingDimension > 1 && obj.timeDelay == 0)
                ERR_MSG = strcat(ERR_MSG, " Value must be greater than", ...
                    " or equal to 1.");
                error(ERR_MSG);
            end
            %}
        end % END set.timeDelay()
        
        function obj = set.threshold(obj, newThreshold)
        %SET.THRESHOLD Set the value of the threshold
        %   Validation: threshold can be any real number greater than or equal 
        %   to 0. Input can be of the form of a single real positive number
        %   (e.g., 0.5) or an interval (e.g., [0.1, 0.25]).
        % -------------------------------------------------------------------- %
            obj.threshold = newThreshold;
        
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
            
            if(obj.threshold(1) < 0 || obj.threshold(2) < 0)
                ERR_MSG = strcat(ERR_MSG, " Value cannot be negative.");
                error(ERR_MSG);
            end
        end % END set.threshold()
        
        function obj = set.normType(obj, newNormType)
        %SET.NORMTYPE Set the norm
        %   Validation: only the norms L1, L2 and L-infinity are accepted. User 
        %   can use these norms as follows:
        %       * as strings: "L1", "L2", and "L3";
        %       * as chars: 'L1', 'L2', and 'L3';
        %       * as numbers: 1, 2, and 3 (respectively to L1, L2, and L3).
        %   Norms are standardized. So, if the user gives "l1", the norm is
        %   changed to 'L1'; in the same way, if the user gives 1, then the norm
        %   is changed to 'L1'.
        %   For L1, user can write: 1, 'l1', 'L1', 'taxicab', "l1", "L1", and
        %   "taxicab"; for L2: 2, 'l2', 'L2', 'euclidean', 'Euclidean', "l2",
        %   "L2", "euclidean", "Euclidean"; for L-infinity: 3, 'l-inf',
        %   'l-infinity', 'L-inf', 'L-infinity', "l-inf", "l-infinity", "L-inf",
        %   "l-infinity", "
        % -------------------------------------------------------------------- %
            obj.normType = newNormType;
            
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
                'l-inf', 'L-inf', 'L-infinity', 'l-infinity', 'Linf', 'linf', ...
                'Chebychev'};
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
                case {'l-infinity', 'L-infinity', 'l-inf', 'L-inf', ...
                        'linf', 'Linf', 'Chebychev'}
                    obj.normType = 'L-infinity';
            end
        end % END set.normType()
        
        function obj = set.data(obj, newData)
            obj.data = newData;
        end % END set.data()
        
    end % END public methods
    
    methods (Access = private)
        %
        % Check data
        %
        
        function chkdata(obj)
        %CHKDATA Check data from time series
        % -------------------------------------------------------------------- %
            if(isempty(obj.data))
                ERR_MSG = "Data from the time series cannot be empty.";
                error(ERR_MSG);
            end
            
            [~, n] = size(obj.data);
            
            if(n == 1)
                % One time series
                RP.chkdatavector(obj.data{1});
            elseif(n == 2)
                % Two times series
                RP.chkdatavector(obj.data{1});
                RP.chkdatavector(obj.data{2});
            else
                ERR_MSG = "Too many time series.";
                error(ERR_MSG);
            end
            
        end % END chkdata()
        
        %
        % Distance matrix
        %
        
        function D = dm(obj)
        %DM Distance matrix
        % -------------------------------------------------------------------- %
            
        
            [~, n] = size(obj.data);
            
            switch n
                case 1
                    % Number of state-space vectors
                    N = numel(obj.data{1}) - obj.timeDelay ...
                        * (obj.embeddingDimension - 1);
                    if(N < 1)
                        ERR_MSG = "Time series is too short.";
                        error(ERR_MSG);
                    end
                case 2
                    % Number of state-space vectors for each time series
                    N1 = numel(obj.data{1}) - obj.timeDelay ...
                        * (obj.embeddingDimension - 1);
                    N2 = numel(obj.data{2}) - obj.timeDelay ...
                        * (obj.embeddingDimension - 1);
                    if(N1 < 1 || N2 < 1)
                        ERR_MSG = "Time series is too short.";
                        error(ERR_MSG);
                    end
            end
            
            % Obtain the self-distance / cross distance matrix
            switch n
                case 1  % If there's only one time series
                    D = sdm(obj);
                case 2  % If there are two time series
                    D = cdm(obj);
            end
        end % END dm()
        
        function D = sdm(obj)
        %SDM Self-distance matrix
        % -------------------------------------------------------------------- %
            % Number of state-space vectors
            N = numel(obj.data{1}) ...
                - obj.timeDelay*(obj.embeddingDimension - 1);

            % State-space vectors (each row of X is a state-space
            % vector)
            X = zeros(N, obj.embeddingDimension);
            for i = 1:obj.embeddingDimension
                X(1:end, i) = obj.data{1}((1:N) + obj.timeDelay*(i - 1));
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

        function D = cdm(obj)
        %CDM Cross-distance matrix
        % -------------------------------------------------------------------- %
            % Number of state-space vectors for each time series
            N1 = numel(obj.data{1}) ...
                - obj.timeDelay*(obj.embeddingDimension - 1);
            N2 = numel(obj.data{2}) ...
                - obj.timeDelay*(obj.embeddingDimension - 1);

            % State-space vectors (each row of X and Y is a state-space vector)
            X = zeros(N1, obj.embeddingDimension);
            Y = zeros(N2, obj.embeddingDimension);
            for i = 1:obj.embeddingDimension
                X(1:end, i) = obj.data{1}((1:N1) + obj.timeDelay*(i - 1));
                Y(1:end, i) = obj.data{2}((1:N2) + obj.timeDelay*(i - 1));
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
    
    methods (Access = private, Static = true)
        %
        % Check data
        %
        
        function chkdatavector(data)
            % Is empty?
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
        end % END chkdatavector()
        
    end % END private, static methods
        
end
