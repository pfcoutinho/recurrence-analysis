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
        dHist   % Histogram of diagonal lines
        vHist   % Histogram of vertical lines
    end
    
    properties (Dependent)
        % Based on the density of recurrence points
        RR                  % Recurrence rate
        
        % Based on diagonal lines
        DET                 % Determinism
        LMAX                % Maximum diagonal line length
        LAVG                % Average diagonal line length
        ENTD                % Entropy of diagonal lines
        DIV                 % Divergence
        RATIO               % Ratio
        
        % Based on vertical lines
        LAM                 % Laminarity
        TT                  % Trapping time
        VMAX                % Maximum vertical line length
        ENTV                % Entropy of vertical lines
        
        % Kolmogorov-Chaitin complexity (compression)
        ZIPD                % ZIP compression of horizontal lines
        ZIPV                % ZIP compression of vertical lines
        TIFF                % TIFF image compression
    end
    
    methods
        %
        % Class constructor function
        %
        
        function obj = RecurrenceAnalysis()
            obj.dHist = hod(RP);
            obj.vHist = hov(RP);
        end
        
    end % END methods
    
    methods (Static)
        %
        % Histograms
        %

        function H = hod(RP)
        %HOD Histogram of diagonal lines
        % -------------------------------------------------------------------- %
            [m, n] = size(RP);
            
            H = zeros(max(m, n), 1);
            
            for i = -(n-1):1:(n-1)
                diagonalLine = [0; find(~diag(RP, i)); numel(diag(RP, i)) + 1];
                lineLengths  = diff(diagonalLine) - 1;
                idx          = (lineLengths > 0);
                H(lineLengths(idx)) = H(lineLengths(idx)) + 1;
            end
        end
        
        function H = hov(RP)
        %HOD Histogram of vertical lines
        % -------------------------------------------------------------------- %
            [m, n] = size(RP);
            
            H = zeros(max(m, n), 1);
            
            for i = 1:m
                verticalLine = [0; find(~RP(i, 1:end)); ...
                    numel(RP(i, 1:end)) + 1];
                lineLengths  = diff(verticalLine) - 1;
                idx          = (lineLengths > 0);
                H(lengths(idx)) = H(lineLengths(idx)) + 1;
            end
        end
        %
        % Recurrence metrics
        %
        
        % Recurrence rate (RR)
        function RR = recurrencerate(RP)
        %RECURRENCERATE Recurrence rate
        % -------------------------------------------------------------------- %
            [m, n] = size(RP);
            % Discard the LOI?
            if(m == n && sum(diag(RP)) == m)
                RR = (nnz(RP) - m)/(m*n);
            else
                RR = nnz(RP)/(m*n);
            end
        end
        
        % Determinism (DET)
        function determinism()
            
        end
        
        function getdeterminism()
        end
        
        % Entropy (ENT)
        function entropyd()
        end
        
        function getentropyd()
        end
        
        function entropyv()
        end
        
        function getentropyv()
        end
        
        % Adaptive threshold
        function adaptivethreshold(rate)
            
        end
    end
end
