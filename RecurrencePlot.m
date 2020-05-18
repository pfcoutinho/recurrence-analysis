classdef RecurrencePlot < DistanceMatrix & ...
                            RecurrenceQuantificationAnalysis & ...
                            handle
%RECURRENCEPLOT Recurrence plot class
%   Recurrence plot (RP) is a tool for non-linear data analysis. It is a binary
%   matrix in which the recurrence points are represented by black or colored
%   points. RPs are a kind of thresholded similarity matrix and reveal interes-
%   ting recurrence patterns of dynamical systems. RP is also knwon as recurren-
%   ce matrix.
%   This class can be used to generate an RP from a single data series or the
%   cross-recurrence plot (CRP) between two data series.
%
% DEPENDECIES
%   DistanceMatrix.m
%
% SYNTAX
%   obj = RecurrencePlot(embeddingDimension, timeDelay, threshold, normType, ...
%           data);
%   obj = RecurrencePlot(embeddingDimension, timeDelay, threshold, normType, ...
%           data1, data2);
%
% INPUT
%   embeddingDimension
%   timeDelay
%   threshold
%   normType
%
% REFERENCES
%   
%
% AUTHOR
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: May 18, 2020
% ============================================================================ %
    %
    % Properties
    %
    properties
        
        threshold                   % threshold (it is estimated if the value of
                                    % the recurrence rate is given)
                                    
        recurrenceRate              % recurrence rate (it is estimated if the 
                                    % the value of the threshold is given)
                                    
        absoluteRecurrenceRateError % absolute recurrence rate error (if the
                                    % recurrenceParameter is the recurrence
                                    % rate, then it is not always possible to
                                    % guarantee that that specific recurrence
                                    % rate will be achieved, so an error term is
                                    % necessary)
        
        RP = []                     % recurrence or cross recurrence plot
                                    
    end %END properties
    
    properties (Access = private)
        
        recurrenceParameter         % keeps track about the recurrence parameter
                                    % being used, if it is the threshold or the
                                    % recurrence rate; recurrenceParameter is
                                    % updated to match the last parameter set by
                                    % user
        
        flagFromInside = false;     % tells if the parameter change comes from 
                                    % within the class (for example, when we are
                                    % updating the value of the threshold after 
                                    % obtaining the RP for a given value of the 
                                    % recurrence rate and vice-versa) or from 
                                    % the outside (for example, when the user 
                                    % set a new value for the threshold or the
                                    % recurrence rate)
    end
    
    %
    % Events
    %
    events (NotifyAccess = private)
        thresholdChangeEvt           % triggered when threshold value changes
        recurrenceRateChangeEvt      % triggered when recurrenceRate value changes
        recurrenceParameterChangeEvt % triggered when recurrenceParameter changes
    end
    
    %
    % Methods
    %
    methods

        function self = RecurrencePlot(embeddingDimension, timeDelay, ...
                            recurrenceParameter, recurrenceParameterValue, ...
                            normType, varargin)
        %RECURRENCEPLOT Recurrence plot class constructor
        % -------------------------------------------------------------------- %
            self = self@DistanceMatrix(embeddingDimension, timeDelay, ...
                        normType, varargin{1:end});
            
            self.recurrenceParameter = string(recurrenceParameter);
            switch self.recurrenceParameter
                case "threshold"
                    self.threshold = recurrenceParameterValue;
                    self = thresholdrp(self);
                case "recurrence rate"
                    self.recurrenceRate = recurrenceParameterValue;
                    self = recurrenceraterp(self);
                otherwise
                    error("Unrecognized recurrence parameter: %s.", ...
                            self.recurrenceParameter);
            end
            
            % these are the listeners to the event that is triggered when
            % threshold or recurrence parameter is set or changes
            addlistener(self, 'thresholdChangeEvt', @thresholdchangeevthandler);
            addlistener(self, 'recurrenceRateChangeEvt', ...
                    @recurrenceratechangeevthandler);
        end %END RecurrencePlot() 
        

        function plot(self, varargin)
        %PLOT Plot recurrence plot
        %   Plot the recurrence or the cross-recurrence plot. These plots are
        %   binary matrix.
        % -------------------------------------------------------------------- %
            imagesc(self.RP);
            
            if(nargin == 2)
                if(strcmp(varargin{1}, 'color'))
                    colormap([1, 1, 1; rand(1, 3).*[0.9, 0.7, 0.95]])
                end
            else
                colormap([1, 1, 1; 0, 0, 0])
            end
            
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
        end %END plot()
        
        
        function set.threshold(self, value)
        %SET.THRESHOLD Assigns threshold parameter
        % -------------------------------------------------------------------- %
            value = self.chkthreshold(value);
            
            % set threshold value and notify listeners if value changes
            if isempty(self.threshold) || self.threshold ~= value
                self.threshold = value;
                notify(self, 'thresholdChangeEvt');
            end
        end %END set.threshold()
        
        function set.recurrenceRate(self, value)
        %SET.RECURRENCERATE Assigns recurrence rate parameter
        % -------------------------------------------------------------------- %
            self.chkrecurrencerate(value);
            
            % Set threshold value and notify listeners if value changes
            if isempty(self.recurrenceRate) || self.recurrenceRate ~= value
                self.recurrenceRate = value;
                notify(self, 'recurrenceRateChangeEvt');
            end
        end %END set.recurrenceRate()
        
    end %END public methods
    
    
    methods (Access = protected)
        
        function self = thresholdrp(self)
        %THRESHOLDRP Calculates the (cross) recurrence plot applying a threshold
        %
        % -------------------------------------------------------------------- %
            [m, n] = size(self.threshold);
            if(m == 1 && n == 1)
                % recurrence plot / cross-recurrence plot
                self.RP = (self.M <= self.threshold);
            else
                % corridor recurrence plot / cross recurrence plot
                self.RP = and(self.M >= self.threshold(1), ...
                        self.M <= self.threshold(2));
            end
            
            % updates the value of the rrecurrenceRate property (without
            % recalculating the RP); also, considers that the
            % absoluteRecurrenceRateError is equal to 0 (because recurrenceRate
            % is computed directly from the RP)
            self.flagFromInside = true;
            self.recurrenceRate = self.rr();
            self.absoluteRecurrenceRateError = 0;
            self.flagFromInside = false;
        end %END rpusingthreshold()
        
        function self = recurrenceraterp(self)
        %RECURRENCERATERP This function employes recurrence rate to calculate
        %the (cross) recurrence plot
        %   It is not always possible to find the threshold that corresponds to
        %   a specific value of the recurrence rate for any given data series.
        %   This is related to structures which are present in recurrence
        %   plots.
        % -------------------------------------------------------------------- %
            [~, n] = size(self.data);
            
            % range of search
            minThreshold = 0;       % default value
            values = [0, 0];
            for i = 1:n
                values(i) = max(self.data{i});
            end
            maxThreshold = max(values);
            
            % binary search
            newThreshold = (minThreshold + maxThreshold)/2;
            while true
                self.RP = (self.M <= newThreshold);
                currentRecurrenceRate = self.rr();
                
                if abs(maxThreshold - minThreshold) < 1e-06 || ...
                        abs(currentRecurrenceRate - self.recurrenceRate) < 1e-06
                    break;
                elseif currentRecurrenceRate > self.recurrenceRate
                    maxThreshold = newThreshold;
                    newThreshold = (minThreshold + maxThreshold)/2;
                elseif currentRecurrenceRate < self.recurrenceRate
                    minThreshold = newThreshold;
                    newThreshold = (minThreshold + maxThreshold)/2;
                end
            end
            
            self.absoluteRecurrenceRateError = ...
                abs(self.recurrenceRate - currentRecurrenceRate);
            
            % updates the value of the threshold property (without recalculating
            % the RP)
            self.flagFromInside = true;
            self.threshold = newThreshold;
            self.flagFromInside = false;
        end %END recurrenceraterp()
        
        
        function self = parameterchangeevthandler(self, ~)
        % -------------------------------------------------------------------- %
            self = parameterchangeevthandler@DistanceMatrix(self);            
            if ~isempty(self.RP)
                switch self.recurrenceParameter
                    case "threshold"
                        self = thresholdrp(self);
                    case "recurrence rate"
                        self = recurrenceraterp(self);
                end
            end
        end
        
        function self = thresholdchangeevthandler(self, ~)
        %
        %   When flagFromInside is UP (i.e., it is equal to 1), the RP isn't
        %   recalculated and recurrenceParameter isn't updated
        % -------------------------------------------------------------------- %
            if ~self.flagFromInside
                if ~isempty(self.RP)
                    % computes a new RP
                    self = thresholdrp(self);
                end
            
                if strcmp(self.recurrenceParameter, "recurrence rate")
                    self.recurrenceParameter = "threshold";
                end
            end
        end
        
        function self = recurrenceratechangeevthandler(self, ~)
        % -------------------------------------------------------------------- %
            if ~self.flagFromInside
                if ~isempty(self.RP)
                    % computes a new RP
                    self = recurrenceraterp(self);
                end

                if strcmp(self.recurrenceParameter, "threshold")
                    self.recurrenceParameter = "recurrence rate";
                end
            end
        end
        
        
        function value = chkthreshold(self, value)
        %CHKTHRESHOLD Validates the threshold parameter
        %
        %   Threshold must be:
        %       * non-empty,
        %       * numeric,
        %       * real valued,
        %       * greater than or equal to 0, and
        %       * a scalar or an interval (e.g., [0.2, 0.4]).
        % -------------------------------------------------------------------- %
            parameterName = "threshold";
            self.verifyifempty(value, parameterName);
            self.verifyifnumeric(value, parameterName);
            self.verifyifreal(value, parameterName);

            % check if value is a scalar or an array; also, if it is an array,
            % check its size and sort the elements
            if ~isscalar(value)
                if ~isvector(value)
                    ERR_MSG = strcat(DEFAULT_ERR_MSG, " It must be a ", ...
                    "scalar or a vector with exactly two elements.");
                    error(ERR_MSG);
                else
                    % check vector sizes
                    [m, ~] = size(value);
                    if m == 2
                        value = value';
                    end
                    
                    % sort elements in ascending order
                    value = sort(value, 'ascend');
                end
            end
        end %END chkthreshold()
        
        function value = chkrecurrencerate(self, value)
        %CHKRECURRENCERATE Validates the recurrence rate parameter
        %
        %   Recurrence rate must be:
        %       * non-empty,
        %       * numeric
        %       * real valued,
        %       * a scalar in the interval [0, 1].
        % -------------------------------------------------------------------- %
            parameterName = "recurrence rate";
            self.verifyifempty(value, parameterName);
            self.verifyifnumeric(value, parameterName);
            self.verifyifscalar(value, parameterName);
            self.verifyifreal(value, parameterName);

            if (value < 0) || (value > 1)
                error(strcat("Invalid parameter: ", parameterName, ...
                    ". Value must be in the interval [0, 1]."));
            end
        end %END chkrecurrencerate()
        
    end %END protected methods
    
    
    methods (Access = protected)
        
        function recurrenceRate = rr(self)
            recurrenceRate = rr@RecurrenceQuantificationAnalysis(self);
        end
        
    end
    
    methods (Access = protected, Static = true)
        
        function verifyifempty(value, parameterName)
            verifyifempty@DistanceMatrix(value, parameterName);
        end %END verifyifempty()
        
        function verifyifnumeric(value, parameterName)
            verifyifnumeric@DistanceMatrix(value, parameterName);
        end %END verifyifnumeric()
        
        function verifyifreal(value, parameterName)
            verifyifreal@DistanceMatrix(value, parameterName);
        end %END verifyifreal()
        
    end %END protected, static methods
end
