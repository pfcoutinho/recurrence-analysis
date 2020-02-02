classdef RecurrencePlot < Recurrence & RQA
%RECURRENCEPLOT Recurrence plot
%   This class can be used to generate the recurrence plot of a time series or
%   the cross-recurrence plot between two time series (in this case, it is not
%   necessary for the time series to have the same length).
%
%   RecurrencePlot class extends Recurrence class and inherits its properties:
%   embeddingDimension, timeDelay, and normType.
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
% PROPERTIES
%   RP
%       Recurrence or cross-recurrence plot (binary matrix)
%
%   See Recurrence class for a detailed description of the other properties.
%
% METHODS
%   plotr()
%       Plot the recurrence or the cross-recurrence plot
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 30, 2020
% ============================================================================ %

    properties
        RP              % Recurrence plot
    end % END properties

    methods
        %
        % Class constructor function
        %

        function obj = RecurrencePlot(embeddingDimension, timeDelay, ...
                threshold, normType, varargin)
            
            % Superclass constructor
            obj = obj@Recurrence(embeddingDimension, timeDelay, threshold, ...
                normType);

            % Obtain the recurrence plot or the cross-recurrence plot
            obj.RP = rp(obj, varargin{1:end});
            
        end % END RecurrencePlot()
        
        %
        % Recurrence plot
        %
        function RP = rp(obj, varargin)
            % Obtain the distance matrix (DM)
            DM = dm(obj, varargin{1:end});
            
            if(obj.threshold(1) == 0)
                RP = (DM <= obj.threshold(2));
            else
                RP = and(DM >= obj.threshold(1), DM <= obj.threshold(2));
            end
        end
        
        %
        % Plot
        %
        
        function plotr(obj)
        %PLOTR Plot recurrence object
        %   In this case, plotr() will plot the recurrence plot or the 
        %   cross-recurrence plot
        % -------------------------------------------------------------------- %
            plotr@Recurrence(obj.RP);
        end % END plotr()
        
        %
        % Set
        %
        
        function obj = setthreshold(obj, newThreshold)
            obj = setthreshold@Recurrence(obj, newThreshold);
        end % END setthreshold()
        
    end % END public methods
end

