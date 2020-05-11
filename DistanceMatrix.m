classdef DistanceMatrix < handle
%DISTANCEMATRIX
%   
%
% DEPENDENCIES
%   None.
%
% SYNTAX
%   obj = DM(embeddingDimension, timeDelay, normType, data);
%   obj = DM(embeddingDimension, timeDelay, normType, data1. data2);
%
% INPUT
%   embeddingDimension       
%   timeDelay
%   normType
%   data
%
% AUTHOR
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% SEE ALSO
%   RecurrencePlot.m
%
% Last update: May 11, 2020
% ============================================================================ %
    %
    % Properties
    %
    properties
        % recurrence parameters
        embeddingDimension          % embedding dimension
        timeDelay                   % time delay
        normType                    % norm (L1, L2, or L-infinity)
        data                        % data from time series
        M = []                      % distance matrix
    end
    
    %
    % Events
    %
    events
        parameterChangeEvt          % event to handle a parameter change
    end
    
    %
    % Methods
    %
    methods
        
        function self = DistanceMatrix(embeddingDimension, timeDelay, ...
                            normType, varargin)
        %DISTANCEMATRIX Distance matrix
        %   This is a function for data series analysis. It calculates the
        %   self-distance matrix respective to a single data series or the 
        %   cross-distance matrix relative to two data series.
        % -------------------------------------------------------------------- %
            switch nargin
                case {4, 5}
                    self.embeddingDimension = embeddingDimension;
                    self.timeDelay          = timeDelay;
                    self.normType           = normType;
                    self.data               = varargin(1:end);
                otherwise
                    error("");
            end
            
            % obtain the distance matrix
            self.M = dm(self);
            
            % add listeners to parameter change events
            addlistener(self, 'parameterChangeEvt', @evthandle);
        end % END DistanceMatrix() 
        
        
        function plot(self)
        %PLOT Plots the self-distance matrix or the cross-distance matrix
        % -------------------------------------------------------------------- %
            imagesc(self.M)
            
            colormap('parula')
            colorbar
            
            [m, n] = size(self.M);
            if(m == n)  % symmetric matrix
                axis square
            end

            % axis
            [m, n] = size(self.M);
            ax = [1, round(n/2), n];
            ay = [1, round(m/2), m];
            set(gca, 'XTick', ax)
            set(gca, 'YTick', ay)
        end % END plot()
        
        
        function set.embeddingDimension(self, value)
        %SET.EMBEDDINGDIMENSION Embedding dimension
        %   Embedding dimension is a recurrence parameter. It is used to
        %   generate the state-space vector and indicates how many samples each
        %   state-space vector contains.
        %
        %   VALIDATION Embedding dimension must be
        %       * non-empty,
        %       * numeric,
        %       * a scalar,
        %       * an integer,
        %       * greater than or equal to 1.
        % -------------------------------------------------------------------- %
            self.chkembeddingdimension(value);
            
            % set embeddingDimension and notify listeners if value changes
            if(isempty(self.embeddingDimension) ...
                    || self.embeddingDimension ~= value)
                self.embeddingDimension = value;
                notify(self, 'parameterChangeEvt');
            end
        end % END set.embeddingDimension()
        
        
        function set.timeDelay(self, value)
        %SET.TIMEDELAY Time delay
        %   Time delay is a recurrence parameter. It is used along with the 
        %   embedding dimension to generate the state-space vector.
        %
        %   VALIDATION Time delay must be:
        %       * non-empty,
        %       * numeric,
        %       * a scalar,
        %       * an integer,
        %       * greater than or equal to 1 if embedding dimension is greater
        %           than 1; if embedding dimension is equal to 1, the value of
        %           time delay is irrelevant
        % -------------------------------------------------------------------- %
            self.chktimedelay(value);
            
            % Set timeDelay and notify listeners if value changes
            if(isempty(self.timeDelay) || self.timeDelay ~= value)
                self.timeDelay = value;
                notify(self, 'parameterChangeEvt');
            end
        end % END set.timedelay()
        
        
        function set.normType(self, value)
        %SET.NORMTYPE Norm
        %
        %   VALIDATION Only the norms L1, L2 and L-infinity are accepted. User 
        %   can use these norms as follows:
        %       * as strings: "L1", "L2", and "L3";
        %       * as chars: 'L1', 'L2', and 'L3';
        %       * as numbers: 1, 2, and 3 (respectively to L1, L2, and L3).
        % -------------------------------------------------------------------- %
            value = self.chknorm(value);
            
            % set normType and notify listeners if value changes
            if(isempty(self.normType) || ~strcmp(self.normType, value))
                self.normType = value;
                notify(self, 'parameterChangeEvt');
            end
        end % END set.normtype()
        
        
        function set.data(self, data)
        %SET.DATA Data from time series
        %
        %   VALIDATION Data from time series must be:
        %       * non-empty
        %       * numeric
        % -------------------------------------------------------------------- %
            self.chkdata(data);
            
            % set data and notify listeners if value changes
            if(isempty(self.data) || ~isequal(self.data, data))
                self.data = data;
                notify(self, 'parameterChangeEvt');
            end
        end % END set.data()
        
    end % END public methods
    
    
    methods (Access = protected)
        
        function D = dm(self)
        %DM Distance matrix
        %   Obtains the self-distance or the cross distance matrix.
        % -------------------------------------------------------------------- %
            [~, n] = size(self.data);
            
            % first, we obtain the state-space vectors for each time series (N
            % being the number of state-space vectors)
            [S, N] = generatessvectors(self);
            
            if n == 1       % only one data series
                % replicate state-space vectors if number of state space vectors
                % is less than 1e5. In this case, computations are faster.
                if N < 1e5
                    X{1} = repmat(S{1}, N, 1);
                    X{2} = reshape(repmat(S{1}, 1, N)', ...
                                        self.embeddingDimension, ...
                                        N*N)';
                    D = calcdistances(self, X{1}, X{2});
                    D = reshape(D, N, N);
                else
                    D = zeros(N, N);
                    for i = 1:N
                        D(i, :) = calcdistances(self, S{1}(i, :), S{1});
                    end
                end
            elseif n == 2   % two data series
                if N(1) < 1e5 && N(2) < 1e5
                    X{1} = repmat(S{1}, N(2), 1);
                    X{2} = reshape(repmat(S{2}, 1, N(1))', ...
                                    self.embeddingDimension, ...
                                    N(1)*N(2))';
                    D = calcdistances(self, X{1}, X{2});
                    D = reshape(D, N(1), N(2));
                else
                    D = zeros(N(1), N(2));
                    if N(1) <= N(2)
                        for i = 1:N(1)
                            D(i, :) = calcdistances(self, S{1}(i, :), S{2});
                        end
                    else
                        for i = 1:N(2)
                            D(:, i) = calcdistances(self, S{2}(i, :), S{1});
                        end
                    end
                end
            end
            
            %
            % auxiliary subfunctions
            %
            function [S, N] = generatessvectors(self)
            %GENERATESSVECTORS Generate state-space vectors
            % ---------------------------------------------------------------- %
                S = cell(1, n);     % state-space vectors
                N = zeros(1, n);    % number of state-space vectors to
                                    % generate for each time series    
                for j = 1:n
                    N(j) = numel(self.data{j}) ...
                            - self.timeDelay*(self.embeddingDimension-1);
                    
                    S{j} = zeros(N(j), self.embeddingDimension, 'single');
                    for k = 1:self.embeddingDimension
                        S{j}(1:end, k) = self.data{j}((1:N(j)) ...
                                    + self.timeDelay*(k-1));
                    end
                end
            end %END generatespvectors()
                
            function D = calcdistances(self, X, Y)
            %CALCDISTANCES
            % ---------------------------------------------------------------- %
                switch self.normType
                    case {'L1', 'l1', 'l-1 norm', 'taxicab'}
                        D = sum(abs(X - Y), 2);
                    case {'L2', 'l2', 'l-2 norm', 'Euclidean', 'euclidean'}
                        D = sqrt(sum((X - Y).^2, 2));
                    case {'L-infinity', 'l-infinity', 'l-infinity norm', 'L-inf', 'l-inf', ...
                            'Maximum', 'maximum'}
                        D = max(abs(X - Y), [], 2);
                end
            end %END calcdistances()
        end % END dm()
        

        function self = evthandle(obj, ~)
        %EVTHANDLE Function to handle the change of a parameter value
        %   It calculates the distance matrix again.
        % -------------------------------------------------------------------- %
            if(~isempty(obj.M))
                self.M = dm(obj);
            end
        end % END evthandle()
        
    end % END protected methods
    
    
    methods (Access = private, Static = true)
        
        function chkembeddingdimension(value)
        %CHKEMBEDDINGDIMENSION 
        % -------------------------------------------------------------------- %
            DEFAULT_ERR_MSG = "Invalid parameter: embedding dimension.";
            
            % check if value is empty
            if(isempty(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value cannot be empty.");
                error(ERR_MSG);
            end
            
            % check if value is numeric
            if(~isnumeric(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value must be numeric.");
                error(ERR_MSG);
            end
            
            % check if value is scalar
            if(~isscalar(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value must be a scalar.");
                error(ERR_MSG);
            end
            
            % check if value is an integer
            if(mod(value, 1) ~= 0)
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value must be an integer.");
                error(ERR_MSG);
            end
            
            % check if value is greater than or equal to 1
            if(value < 1)
               ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value must be greater ", ...
                   "than or equal to 1.");
                error(ERR_MSG);
            end
        end
        
        
        function chktimedelay(value)
        %CHKTIMEDELAY
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
        end
        
        
        function value = chknorm(value)
        %CHKNORM Validate the norm
        %   We will check the norm being used and, if possible, we will
        %   standardize it.
        % -------------------------------------------------------------------- %
            DEFAULT_ERR_MSG = "Invalid parameter: norm.";
            
            % Check if value is empty
            if(isempty(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Norm cannot be empty.");
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
                        value = 'l-1 norm';
                    case 2
                        value = 'l-2 norm';
                    case 3
                        value = 'l-infinity norm';
                    otherwise
                        ERR_MSG = strcat(DEFAULT_ERR_MSG, " If a number, ", ...
                            "value must be 1, 2, or 3, for 'l1', l2, ", ...
                            "or 'l-infinity', respectively.");
                        error(ERR_MSG);
                end
            end
            
            % check if value is char
            if(~ischar(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Unrecognizable norm.");
                error(ERR_MSG);
            end
            
            % check if value is a valid norm
            VALID_NORMS = {'L1','l1', 'L-1', 'l-1', 'l1 norm', 'l-1 norm', ...
                'taxicab', 'L2', 'l2', 'L-2', 'l-2', 'l2 norm', 'l-2 norm', ...
                'L2 norm', 'L-2 norm', 'euclidean', 'Euclidean', ...
                'L-inf', 'l-inf', 'l-inf norm', 'L-inf norm', 'L-infinity', ...
                'l-infinity', 'linf', 'Linf', 'Linfinity', 'linfinity', ....
                'Chebychev'};
            if(~ismember(value, VALID_NORMS))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Choose among 'l-1', ", ...
                    "'l-2', and 'l-infinity' norms.");
                error(ERR_MSG);
            end
            
            % standardize the norm
            switch value
                case {'L1','l1', 'L-1', 'l-1', 'l1 norm', 'taxicab'}
                    value = 'l-1 norm';
                case {'L2', 'l2', 'L-2', 'l-2', 'l2 norm', 'l-2 norm', ...
                        'L2 norm', 'L-2 norm', 'euclidean', 'Euclidean'}
                    value = 'l-2 norm';
                case {'L-inf', 'l-inf', 'l-inf norm', 'L-inf norm', ...
                        'L-infinity', 'l-infinity', 'linf', 'Linf', ....
                        'Linfinity', 'linfinity', 'Chebychev'}
                    value = 'l-infinity norm';
            end
        end
        
        
        function data = chkdata(data)
        %CHKDATA Check data from time series
        % -------------------------------------------------------------------- %
            DEFAULT_ERR_MSG = "Invalid time series.";

            % check if data is empty
            if(isempty(data))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Time series cannot ", ...
                    "be empty.");
                error(ERR_MSG);
            end
            
            % data must be a cell
            if(~iscell(data))
                data = {data};
            end
                
            [~, n] = size(data);
            
            if n > 2
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Too many data series.");
                error(ERR_MSG);
            end
            
            for i = 1:n
                % check if data is a vector
                if(~isvector(data{i}))
                    ERR_MSG = strcat(DEFAULT_ERR_MSG, "Time series must be a ", ...
                        "vector.");
                    error(ERR_MSG);
                end
                
                % check if it is numeric or logical
                if(~isnumeric(data{i}) && ~islogical(data{i}))
                    ERR_MSG = "Time series must be numeric or logical.";
                    error(ERR_MSG);
                end 
            end
        end % END chkdata()
        
    end % END private, static methods
    
end
