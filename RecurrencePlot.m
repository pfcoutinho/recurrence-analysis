classdef RecurrencePlot < Recurrence
%RECURRENCEPLOT Recurrence plot
%   This class can be used to generate the recurrence plot of a time series or
%   the cross-recurrence plot between two time series (in this case, it is not
%   necessary for the time series to have the same length).
%
%   RecurrencePlot class extends Recurrence class and inherits its properties:
%   embeddingDimension, timeDelay, and normType.
%
% PROPERTIES
%   RP
%       Recurrence or cross-recurrence plot (binary matrix)
%
%   See Recurrence class for a description of the other properties.
%
% METHODS
%   plotr()
%       Plot the recurrence or the cross-recurrence plot
%
% SYNTAX
%   To obtain the recurrence plot:
%   RPobj = RecurrencePlot(embeddingDimension, timeDelay, threshold, normType,
%       timeSeries)
%
%   To obtain the cross-recurrence plot:
%   RPobj = RecurrencePlot(embeddingDimension, timeDelay, threshold, normType,
%       timeSeries1, timeSeries2)
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 30, 2020
% ============================================================================ %

    properties
        RP              % Recurrence plot
    end

    methods
        %
        % Class constructor function
        %

        function obj = RecurrencePlot(embeddingDimension, timeDelay, ...
                threshold, normType, varargin)
            
            % Superclass constructor
            obj = obj@Recurrence(embeddingDimension, timeDelay, threshold, ...
                normType);

            % Obtain the distance matrix (DM)
            obj.RP = dm(obj, varargin{1:end});
            
            % Apply the threshold
            recurrenceType = checkrecurrencetype(obj);
            
            switch recurrenceType
                case 'normal'
                    obj.RP = (obj.RP <= obj.threshold);
                case 'corridor'
                    obj.RP = and(obj.RP >= obj.threshold(1), ...
                        obj.RP <= obj.threshold(2));
            end
            
            % Make it sparse? Yes, if the quantity of recurrence points is
            % less than 20% of the total points
            [m, n] = size(obj.RP);
            if(sum(obj.RP(1:end)/(m*n)) < 0.2)
                obj.RP = sparse(obj.RP);
            end
            
        end % END RecurrencePlot()
        
        %
        % Plot
        %
        
        function plotr(obj)
        %PLOTR
        % -------------------------------------------------------------------- %
            plotr@Recurrence(obj.RP);
        end
        
        %
        % Set
        %
        
        function obj = setthreshold(obj, newThreshold)
            obj = setthreshold@Recurrence(obj, newThreshold);
        end
        
    end % END public methods
    
    methods (Access = protected)
        
        function recurrenceType = checkrecurrencetype(obj)
        %CHECKRECURRENCETYPE
        % -------------------------------------------------------------------- %
            [m, n] = size(obj.threshold);
            
            if(m == 1 && n == 1)
                recurrenceType = 'normal';
            else
                % Corridor thresholded version of the recurrence plot
                recurrenceType = 'corridor';
            end
        end % END checkrecurrencetype()
        
    end % END protected methods
end

