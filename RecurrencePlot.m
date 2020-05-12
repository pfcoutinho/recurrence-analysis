classdef RecurrencePlot < DistanceMatrix & handle
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
% Last update: May 11, 2020
% ============================================================================ %
    %
    % Properties
    %
    properties
        threshold           % threshold
        
        RP = []             % recurrence plot / cross recurrence plot
    end %END properties
    
    %
    % Events
    %
    events
        thresholdChangeEvt
    end
    
    %
    % Methods
    %
    methods

        function self = RecurrencePlot(embeddingDimension, timeDelay, ...
                            threshold, normType, varargin)
        %RECURRENCEPLOT recurrence plot
        % -------------------------------------------------------------------- %
            self = self@DistanceMatrix(embeddingDimension, timeDelay, ...
                        normType, varargin{1:end});
            
            % set threshold
            self.threshold = threshold;
            
            % calculates the recurrence matrix
            self.RP = rp(self);

            % add listeners to the event that is triggered when threshold is
            % changed
            addlistener(self, 'thresholdChangeEvt', @thresholdevthandler);
        end %END RecurrencePlot() 
        

        function plot(self, varargin)
        %PLOT Plot recurrence plot
        %   Plot the recurrence or cross-recurrence plot
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
        %SET.THRESHOLD Threshold
        %
        %   VALIDATION Threshold must be:
        %       * non-empty,
        %       * real valued,
        %       * greater than or equal to 0, and
        %       * a scalar or an interval (e.g., [0.2, 0.4]).
        % -------------------------------------------------------------------- %
            value = self.chkthreshold(value);
            
            % Set threshold value and notify listeners if value changes
            if(isempty(self.threshold) || self.threshold ~= value)
                self.threshold = value;
                notify(self, 'thresholdChangeEvt');
            end
        end %END set.threshold()
        
    end %END public methods
    
    
    methods (Access = protected)
        
        function self = evthandler(self, ~)
        % -------------------------------------------------------------------- %
            self = evthandler@DistanceMatrix(self);            
            if ~isempty(self.RP)
                self.RP = rp(self);
            end
        end
        
        function self = thresholdevthandler(self, ~)
        % -------------------------------------------------------------------- %
            if ~isempty(self.RP)
                self.RP = rp(self);
            end
        end
        
        
        function RP = rp(self)
        %RP Recurrence plot
        %   This function applies the threshold (or the corridor threshold) to
        %   the distance matrix.
        % -------------------------------------------------------------------- %
            [m, n] = size(self.threshold);
            
            if(m == 1 && n == 1)
                % recurrence plot / cross-recurrence plot
                RP = (self.M <= self.threshold);
            else
                % corridor recurrence plot / cross recurrence plot
                RP = and(self.M >= self.threshold(1), ...
                        self.M <= self.threshold(2));
            end 
        end %END rp()
        
        
        function value = chkthreshold(self, value)
        %CHKTHRESHOLD Threshold
        %   Validate threshold value before assigning it to a RecurrencePlot
        %   object.
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
        
    end %END protected methods
    
    
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
