classdef RecurrencePlot < DistanceMatrix & handle
%RECURRENCEPLOT Recurrence Plot
%   
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
% AUTHOR
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Mar 10, 2020
% ============================================================================ %
    %
    % Properties
    %
    properties
        % Recurrence parameters
        % Threshold
        threshold double {mustBeNumeric, mustBeReal, mustBeNonnegative}
    end % END properties
    
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
            
            % Recurrence parameters
            self.threshold = threshold;
                        
            % Recurrence matrix
            self.M = rp(self);
        end % END RecurrencePlot() 
        

        function plot(self, varargin)
        %PLOT Plot recurrence plot
        %   Plot the recurrence or cross-recurrence plot
        % -------------------------------------------------------------------- %
            imagesc(self.M)

            if(nargin == 2)
                if(strcmp(varargin{1}, 'color'))
                    colormap([1, 1, 1; rand(1, 3).*[1, 0.7, 1]])
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
        end % END plot()
        

        function set.threshold(self, value)
        %SET.THRESHOLD Set the value of the threshold
        %   Validation: threshold can be any real number greater than or equal 
        %   to 0. Input can be of the form of single scalar (e.g., 0.5) or of
        %   an interval (e.g., [0.1, 0.2]).
        % -------------------------------------------------------------------- %            
            DEFAULT_ERR_MSG = "Invalid parameter: threshold.";
            
            % Check if value is empty
            if(isempty(value))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " Value cannot be empty.");
                error(ERR_MSG);
            end
            
            % Check size of value
            [m, n] = size(value);
                
            if(m < 1 || m > 2 || n < 1 || n > 2 || (m == 2 && n == 2))
                ERR_MSG = strcat(DEFAULT_ERR_MSG, " It must be a scalar ", ....
                    "or a vector with exactly two elements.");
                error(ERR_MSG);
            end
            
            % Check if threshold is already set
            if(isempty(self.threshold))
                flag = false;
            else
                flag = true;
            end

            % Set threshold value and notify listeners if value changes
            if(~flag || self.threshold ~= value)
                self.threshold = value;
                
                if(flag)
                    notify(self, 'parameterChangeEvt');
                end
            end
        end % END set.threshold()
        
    end % END public methods
    
    
    methods (Access = protected)
        
        function self = evthandle(self, ~)
        % -------------------------------------------------------------------- %
            self = evthandle@DistanceMatrix(self);
            self.M = rp(self);
        end
        
    end % END protected methods
    
    
    methods (Access = private)
        
        function M = rp(self)
        % -------------------------------------------------------------------- %
            [m, n] = size(self.threshold);
            
            if(m == 1 && n == 1)
                % Recurrence plot
                M = (self.M <= self.threshold);
            else
                % Corridor recurrence plot
                M = and(self.M >= self.threshold(1), ...
                        self.M <= self.threshold(2));
            end 
        end % END rp()
        
    end % END private methods
    
end
