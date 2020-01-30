classdef RQA
%RQA Recurrence quantification analysis
%   This class performs recurrence quantification analysis (RQA) of recurrence 
%   plots, which are binary matrices (i.e., A(i, j) = {0, 1}).
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 28, 2020
% ============================================================================ %
    
    properties
        diagonalHistogram   % Histogram of diagonal lines
        verticalHistogram   % Histogram of vertical lines
    end
    
    properties (Dependent)
        RR                  % Recurrence rate (RR)
        DET                 % Determinism (DET)
        ENTD                % Entropy of diagonal lines (ENTD)
        ENTV                % Entropy of vertical lines (ENTV)
    end
    
    methods
        function obj = RecurrenceAnalysis()
            
        end
        
        function recurrencerate()
            
        end
        
        function getrecurrencerate()
            
        end
        
        function determinism()
            
        end
        
        function getdeterminism()
        end
        
        function entropyd()
        end
        
        function getentropyd()
        end
        
        function entropyv()
        end
        
        function getentropyv()
        end
        
        function adaptivethreshold(rate)
        end
    end
end
