classdef RecurrencePlot < DistanceMatrix
%RECURRENCEPLOT Recurrence plot
%   This class can be used to generate the recurrence plot of a time series or
%   the cross-recurrence plot between two time series (in this case, it is not
%   necessary for the time series to have the same length).
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 28, 2020
% ============================================================================ %

    properties
        threshold       % Threshold
    end

    methods
        %
        % Class constructor function
        %

        function obj = RecurrencePlot(embeddingDimension, timeDelay, ...
                threshold, normType, varargin)

            % Call superclass
            obj@DistanceMatrix(embeddingDimension, timeDelay, normType, ...
                varargin{1:end});
            
            % Threshold
            obj.threshold = threshold;
            
            if(validateparameter(obj))
            	recurrenceType = checkrecurrencetype(obj);
            end
            
            % Generate the normal / corridor version of the recurrence plot or 
            % the cross-recurrence plot
            switch recurrenceType
                case 'normal'
                    obj.M = (obj.M <= obj.threshold);
                case 'corridor'
                    obj.M = and(obj.M >= obj.threshold(1), ...
                        obj.M <= obj.threshold(2));
            end
            
        end % END RecurrencePlot()
        
        %
        % Validation (threshold parameter)
        %
        
        function returnValue = validateparameter(obj)
        %VALIDATEPARAMETER
        % -------------------------------------------------------------------- %
            returnValue = validatethreshold(obj);
        end
        
        function returnValue = validatethreshold(obj)
        %VALIDATETHRESHOLD
        % -------------------------------------------------------------------- %
            if(isscalar(obj.threshold))

            elseif(isvector(obj.threshold))
                % Corridor thresholded version of the recurrence plot
                obj.threshold = sort(obj.threshold, 'ascend');
            else
                error("Invalid threshold parameter")
            end

            returnValue = true;
        end % END validatethreshold()
        
        
        function recurrenceType = checkrecurrencetype(obj)
        %CHECKRECURRENCETYPE
        % -------------------------------------------------------------------- %
            if(isscalar(obj.threshold))
                recurrenceType = 'normal';
            elseif(isvector(obj.threshold))
                % Corridor thresholded version of the recurrence plot
                recurrenceType = 'corridor';
            end
        end % END checkrecurrencetype()

        %
        % Plot
        %
        
        function plot(obj)
        %PLOT
        % -------------------------------------------------------------------- %
            plot@DistanceMatrix(obj);
            
            % Set colors to B&W
            colormap([1, 1, 1; 0, 0, 0])
        end
        
        %
        % Get and set
        %
        
        function obj = setThreshold(obj, newThreshold)
            obj.threshold = newThreshold;
        end
        
    end % END of methods
end

