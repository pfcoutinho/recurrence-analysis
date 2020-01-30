classdef DistanceMatrix < Recurrence
%DISTANCEMATRIX Distance matrix
%   This class generates the self-distance matrix of a single time series or 
%   the cross-recurrence matrix between two time-series.
%
%   DistanceMatrix class extends Recurrence class and inherits its properties:
%   embeddingDimension, timeDelay, and normType.
%
% PROPERTIES
%   DM
%       Self-distance or cross-distance matrix
%
%   See Recurrence class for a description of the other properties.
%
% METHODS
%   plotr()
%       Plot the self-distance or the cross-distance matrix
%
% SYNTAX
%   To obtain the self-distance matrix:
%   DMobj = DistanceMatrix(embeddingDimension, timeDelay, normType, timeSeries)
%
%   To obtain the cross-distance matrix:
%   DMobj = DistanceMatrix(embeddingDimension, timeDelay, normType, timeSeries1,
%       timeSeries2)
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 30, 2020
% ============================================================================ %

    properties
        DM              % Self-distance or cross-distance matrix
    end

    methods
        %
        % Class constructor
        %
        
        function obj = DistanceMatrix(embeddingDimension, timeDelay, ...
                normType, varargin)
        
            % Superclass constructor
            obj = obj@Recurrence(embeddingDimension, timeDelay, normType);
            
            % Obtain the distance matrix (DM)
            obj.DM = dm(obj, varargin{1:end});
            
        end % END DistanceMatrix()
        
        %
        % Plot
        %
        
        function plotr(obj)
        %PLOTR
        % -------------------------------------------------------------------- %
            plotr@Recurrence(obj.DM);
        end
        
        %
        % Sets
        % 
        
        function obj = setembeddingdimension(obj, newEmbeddingDimension)
        %SETEMBEDDINGDIMENSION
        % -------------------------------------------------------------------- %
            obj.embeddingDimension = newEmbeddingDimension;
        end % END setembeddingdimension()
        
        function obj = settimedelay(obj, newTimeDelay)
        %SETTIMEDELAY
        % -------------------------------------------------------------------- %
            obj = settimedelay@Recurrence(obj, newTimeDelay);
        end % END settimedelay()
        
        function obj = setnormtype(obj, newNormType)
        %SETNORMTYPE
        % -------------------------------------------------------------------- %
            obj.normType = newNormType;
        end % END setnormtype()
        
    end % END public methods
    
    methods (Access = protected)
        %
        % Data validation
        %
        
        function checkdata()
        %CHECKDATA
        % -------------------------------------------------------------------- %
        
        end % END checkdata()
        
        
        
    end % END protected methods
    
end

