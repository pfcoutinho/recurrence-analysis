classdef DistanceMatrix < handle
%DM Distance matrix
%
% SYNTAX
%   obj = DM(embeddingDimension, timeDelay, normType, data);
%   obj = DM(embeddingDimension, timeDelay, normType, data1, data2);
%
% PROPERTIES
%   embeddingDimension
%       Embedding dimension (positive integer greater than or equal to 1)
%
%   timeDelay
%       Time delay (positive integer greater than or equal to 1)
%
%
%   normType
%       Norm used to compute the distance between state-space vectors
%
%   M
%       The recurrence matrix, which is the recurrence plot or the
%       cross-recurrence plot
%
% METHODS
%   DM()
%       Class constructor
%
%   plot()
%       Plot distance matrix
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Mar 4, 2020
% ============================================================================ %
    %
    % Properties
    %
    properties
        % Recurrence parameters
        % Embedding dimension
        embeddingDimension double {mustBeNumeric, mustBeInteger, ...
                                   mustBePositive}
        % Time delay
        timeDelay double {mustBeNumeric, mustBeInteger, ...
                          mustBeNonnegative}

        % Norm (L1, L2, or L-infinity)
        normType
        
        % Data from time series
        timeSeries
        
        % Distance matrix
        M
    end
    
    %
    % Events
    %
    events
        parameterChangeEvt
    end
    
    %
    % Methods
    %
    methods
        %
        % Class constructor
        %
        function self = DistanceMatrix(embeddingDimension, timeDelay, ...
                            normType, varargin)
        %DISTANCEMATRIX Distance matrix
        % -------------------------------------------------------------------- %
            switch nargin
                case {4, 5}
                    self.embeddingDimension = embeddingDimension;
                    self.timeDelay          = timeDelay;
                    self.normType           = normType;
                    self.timeSeries         = varargin(1:end);
                    
                    % Obtain the distance matrix
                    self.M                  = dm(self);
            end
            
            % Add listeners to parameter change event
            addlistener(self, 'parameterChangeEvt', @evthandle);
        end % END DistanceMatrix() 
        
        %
        % Plot
        %
        function plot(self)
        %PLOT Plot distance matrix
        %   Plots the self-distance or the cross distance matrix
        % -------------------------------------------------------------------- %
            imagesc(self.M)
            
            colormap('parula')
            colorbar
            
            [m, n] = size(self.M);
            if(m == n)  % symmetric matrix
                axis square
            end

            % Axis
            [m, n] = size(self.M);
            ax = [1, round(n/2), n];
            ay = [1, round(m/2), m];
            set(gca, 'XTick', ax)
            set(gca, 'YTick', ay)
        end % END plot()
        
        %
        % Sets
        %
        function set.embeddingDimension(self, value)
        %SET.EMBEDDINGDIMENSION Set the value of the embedding dimension
        %   Validation: embedding dimension must be a integer greater than or 
        %   equal to 1.
        % -------------------------------------------------------------------- %
            DEFAULT_ERR_MSG = "Invalid parameter: embedding dimension.";
        
            % Check if value is empty
            if(isempty(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value cannot be empty.");
                error(ERR_MSG);
            end
            
            % Check if value is scalar
            if(~isscalar(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value must be a scalar.");
                error(ERR_MSG);
            end
            
            % Check if embeddingDimension is already set
            if(isempty(self.embeddingDimension))
                flag = false;
            else
                flag = true;
            end
            
            % Set embeddingDimension and notify listeners if value changes
            if(~flag || self.embeddingDimension ~= value)
                self.embeddingDimension = value;

                if(flag)
                    notify(self, 'parameterChangeEvt');
                end
            end
        end % END set.embeddingDimension()
        
        function set.timeDelay(self, value)
        %SET.TIMEDELAY Set the value of the time delay
        %   Validation: the value of the time delay must be greater than or
        %   equal to 1.
        % -------------------------------------------------------------------- %
            DEFAULT_ERR_MSG = "Invalid parameter: time delay.";
        
            % Check if value is empty
            if(isempty(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value cannot be empty.");
                error(ERR_MSG);
            end
            
            % Check if value is scalar
            if(~isscalar(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value must be a scalar.");
                error(ERR_MSG);
            end
            
            % Check if timeDelay is already set
            if(isempty(self.timeDelay))
                flag = false;
            else
                flag = true;
            end
            
            % Set timeDelay and notify listeners if value changes
            if(~flag || self.timeDelay ~= value)
                self.timeDelay = value;
                             
                if(flag)
                    notify(self, 'parameterChangeEvt');
                end
            end
        end % END set.timedelay()
        
        function set.normType(self, value)
        %SET.NORMTYPE Set the norm
        %   Validation: only the norms L1, L2 and L-infinity are accepted. User 
        %   can use these norms as follows:
        %       * as strings: "L1", "L2", and "L3";
        %       * as chars: 'L1', 'L2', and 'L3';
        %       * as numbers: 1, 2, and 3 (respectively to L1, L2, and L3).
        % -------------------------------------------------------------------- %
            DEFAULT_ERR_MSG = "Invalid parameter: norm.";
            
            % Check if value is empty
            if(isempty(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value cannot be empty.");
                error(ERR_MSG);
            end
            
            % Check if value is a string
            if(isstring(value))
                value = char(value);
            end
            
            % Check if value is numeric
            if(isnumeric(value))
                switch(value)
                    case 1
                        value = 'l1';
                    case 2
                        value = 'l2';
                    case 3
                        value = 'l-infinity';
                    otherwise
                        ERR_MSG = strcat(DEFAULT_ERR_MSG, " If a number, ", ...
                            "value must be 1, 2, or 3, for 'l1', l2, ", ...
                            "or 'l-infinity', respectively.");
                        error(ERR_MSG);
                end
            end
            
            % Check if value is char
            if(~ischar(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Unrecognizable norm.");
                error(ERR_MSG);
            end
            
            % Check if value is a valid norm
            VALID_NORMS = {'L1','l1', 'L-1', 'l-1', 'taxicab', ...
                           'L2', 'l2', 'L-2', 'l-2', 'euclidean', 'Euclidean', ...
                           'L-inf', 'l-inf', 'L-infinity', 'l-infinity', ...
                           'linf', 'Linf', 'Linfinity', 'linfinity', ....
                           'Chebychev'};
            if(~ismember(value, VALID_NORMS))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Choose among 'l1', ", ...
                    "'l2', and 'l-infinity'.");
                error(ERR_MSG);
            end
            
            % Check if normType is already set
            if(isempty(self.normType))
                flag = false;
            else
                flag = true;
            end
            
            % Set normType and notify listeners if value changes
            if(~flag || ~strcmp(self.normType, value))
                self.normType = value;

                if(flag)
                    notify(self, 'parameterChangeEvt');
                end
            end
        end % END set.normtype()
        
        function set.timeSeries(self, data)
        %SET.TIMESERIES Set time series
        % -------------------------------------------------------------------- %
            DEFAULT_ERR_MSG = "Invalid time series.";
        
            % Check if data is empty
            if(isempty(data))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Time series cannot ", ...
                    "be empty.");
                error(ERR_MSG);
            end
            
            % Check if data is a cell
            if(~iscell(data))
                data = {data};
            end
            
            % Check size of data
            [~, n] = size(data);
            
            if(n == 1)
                % One time series
                self.chkdata(data{1}, DEFAULT_ERR_MSG);
                
            elseif(n == 2)
                % Two times series
                self.chkdata(data{1}, DEFAULT_ERR_MSG);
                self.chkdata(data{2}, DEFAULT_ERR_MSG);
            elseif(n > 2)
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Too many time series.");
                error(ERR_MSG);
            end
            
            % Check if timeSeries is already set
            if(isempty(self.timeSeries))
                flag = false;
            else
                flag = true;
            end
            
            % Set timeSeries and notify listeners if value changes
            if(~flag || ~isequal(self.timeSeries, data))
                self.timeSeries = data;
                
                if(flag)
                    notify(self, 'parameterChangeEvt');
                end
            end
        end % END set.timeseries()
    end % END public methods
    
    methods (Access = protected)
        %
        % Distance matrix
        %
        function D = dm(self)
        %DM Distance matrix
        %   Obtains the self-distance or the cross distance matrix.
        % -------------------------------------------------------------------- %
            [~, n] = size(self.timeSeries);
            
            switch n
                case 1
                    % Self-distance matrix (one time series)
                    D = sdm(self);
                case 2
                    % Cross distance matrix (two time series)
                    D = cdm(self);
            end
            
            function D = sdm(self)
            %SDM Self-distance matrix
            % ---------------------------------------------------------------- %
            % Number of state-space vectors
                N = numel(self.timeSeries{1}) ...
                    - self.timeDelay*(self.embeddingDimension - 1);

                % State-space vectors (each row of X is a state-space
                % vector)
                X = zeros(N, self.embeddingDimension);
                for i = 1:self.embeddingDimension
                    X(1:end, i) = self.timeSeries{1}((1:N) ....
                                  + self.timeDelay*(i - 1));
                end

                % Replicate state-space vectors
                U = repmat(X, N, 1);
                V = reshape(repmat(X, 1, N)', self.embeddingDimension, N*N)';

                % Compute the distance between the state-space vectors
                switch self.normType
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

            function D = cdm(self)
            %CDM Cross-distance matrix
            % ---------------------------------------------------------------- %
                % Number of state-space vectors for each time series
                N1 = numel(self.timeSeries{1}) ...
                     - self.timeDelay*(self.embeddingDimension-1);
                N2 = numel(self.timeSeries{2}) ...
                     - self.timeDelay*(self.embeddingDimension-1);

                % State-space vectors (each row of X and Y is a state-space vector)
                X = zeros(N1, self.embeddingDimension);
                Y = zeros(N2, self.embeddingDimension);
                for i = 1:self.embeddingDimension
                    X(1:end, i) = self.timeSeries{1}((1:N1) ...
                                  + self.timeDelay*(i-1));
                    Y(1:end, i) = self.timeSeries{2}((1:N2) ...
                                  + self.timeDelay*(i-1));
                end

                % Replicate state-space vectors
                U = reshape(repmat(X, 1, N2)', self.embeddingDimension, N1*N2)';
                V = repmat(Y, N1, 1);

                % Compute the distance between the state-space vectors
                switch self.normType
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
        end % END dm()
        
        %
        % Event handle
        %
        function obj = evthandle(obj, ~)
        %EVTHANDLE Event handle
        %   This function handles the change of a parameter value. It calculates
        %   the distance matrix again.
        % -------------------------------------------------------------------- %
            obj.M = dm(obj);
        end % END evthandle()
    end % END protected methods
    
    methods (Access = private, Static = true)
        %
        % Check data from time series
        %
        function chkdata(data, DEFAULT_ERR_MSG)
            % Check if data vector is empty
            if(isempty(data))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Data from the time ", ... 
                    "series cannot be empty.");
                error(ERR_MSG);
            end
            
            % Check if data is a vector
            if(~isvector(data))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, "Time series must be a ", ...
                    "vector.");
                error(ERR_MSG);
            end

            % Is numeric?
            if(~isnumeric(data) && ~islogical(data))
                ERR_MSG = "Time series must be numeric or logical.";
                error(ERR_MSG);
            end

        end % END chkdata()
    end % END private, static methods
    
end
