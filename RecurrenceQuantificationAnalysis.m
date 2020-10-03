classdef RecurrenceQuantificationAnalysis < handle
%RECURRENCEQUANTIFICATIONANALYSIS Recurrence Quantification Analysis class
%   This class performs recurrence quantification analysis (RQA) of recurrence 
%   plots, which are binary matrices (i.e., A(i, j) = {0, 1}).
%
%   Measures:
%       Histograms:
%         DHIST         Histogram of diagonal lines
%         VHIST         Histogram of vertical lines
%
%       Based on the density of recurrence points
%         RR            Recurrence rate
%
%       Based on diagonal lines
%         DET           Determinism
%         LMAX          Maximum diagonal line length
%         LAVG          Average diagonal line length
%         ENTD          Entropy of diagonal lines
%         DIV           Divergence
%         RATIO         Ratio
%
%       Based on vertical lines
%         LAM           Laminarity
%         TT            Trapping time
%         VMAX          Maximum vertical line length
%         ENTV          Entropy of vertical lines
%
%       Kolmogorov-Chaitin complexity (compression based measures)
%         ZIPD          ZIP compression of horizontal lines
%         ZIPV          ZIP compression of vertical lines
%         TIFF          TIFF image compression
%
% AUTHOR
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% SEE ALSO
%
% Last update: May 18, 2020
% ============================================================================ %
    %
    % Properties
    %
    properties (Abstract = true)
        threshold
        recurrenceRate
        absoluteRecurrenceRateError
        
        RP
    end

    properties (Access = private)
        % Histograms
        diagonalLinesHistogram  % Histogram of diagonal lines
        verticalLinesHistogram  % Histogram of vertical lines
    end
    
    
    %
    % Methods
    %
    methods (Access = protected)
        function recurrenceRate = rr(self)
        %RR Recurrence rate
        %   Recurrence rate is the density of recurrence points. It indicates
        %   the probability of a state recurring.
        % -------------------------------------------------------------------- %        
            [m, n] = size(self.RP);
            recurrenceRate = nnz(self.RP)/(m*n);
        end % END rr()


        %
        % Histograms
        %
        function H = dlhistogram(self)
        %DHIST Histogram of diagonal lines
        % -------------------------------------------------------------------- %
            if isempty(self.diagonalLinesHistogram)
                H = hod(self);
            else
                H = self.diagonalLinesHistogram;
            end
            
            function H = hod(self)
            % ---------------------------------------------------------------- %
                [m, n] = size(self.RP);

                H = zeros(max(m, n), 1);

                for i = -(n-1):1:(n-1)
                    diagonalLine = [0; find(~diag(self.RP, i)); ...
                                        numel(diag(self.RP, i)) + 1];
                    lineLengths  = diff(diagonalLine) - 1;
                    idx          = (lineLengths > 0);
                    H(lineLengths(idx)) = H(lineLengths(idx)) + 1;
                end
            end % END hod()
        end % END dlhist()
        
        
        function H = vlhist(self)
        %VHIST Histogram of vertical lines
        % -------------------------------------------------------------------- %
            if isempty(self.verticalLinesHistogram)
                H = hov(self);
            else
                H = self.verticalLinesHistogram;
            end
            
            function H = hov(self)
            %HOD Histogram of vertical lines
            % ---------------------------------------------------------------- %
                [m, n] = size(self.RP);

                H = zeros(max(m, n), 1);

                for i = 1:m
                    verticalLine = [0; find(~self.RP(1:end, i)); ...
                                        numel(self.RP(1:end, i)) + 1];
                    lineLengths  = diff(verticalLine) - 1;
                    idx          = (lineLengths > 0);
                    H(lineLengths(idx)) = H(lineLengths(idx)) + 1;
                end
            end % END hov()
        end % END vlhist()
        
        %{
        function value = AVGN(obj)
        %AVGN Average number of neighbors
        % -------------------------------------------------------------------- %
            if(isempty(obj.M))
                error("Cannot calculate the average number of neighbors (AVGN).")
            end
            
            [m, n] = size(obj.M);
            
            if(m == n && sum(diag(obj.M) == m))
                value = nnz(obj.M)/(n);
            end
        end % END AVGN()
        
        %
        % Measures based on diagonal lines
        %
        function value = DET(obj, lmin)
        %DET Determinism
        % -------------------------------------------------------------------- %
            if(isempty(obj.diagonalLinesHistogram))
                notify(obj, 'callHistCalc')
            end
            H = obj.diagonalLinesHistogram;
        
            if(isempty(H))
                ERR_MSG = ['Cannot calculate DET: histogram of diagonal ', ...
                    'lines isn''t available.'];
                error(ERR_MSG);
            end
            
            if(~exist('lmin', 'var'))
                lmin = 2;
            end
            
            N = length(H);
                
            value = (lmin:N)*H(lmin:N)/((1:N)*H);
        end
        
        function value = AVGDL(obj, lmin)
        %AVGL Average diagonal line length
            if(isempty(obj.diagonalLinesHistogram))
                notify(obj, 'callHistCalc')
            end
            H = obj.diagonalLinesHistogram;
            
            if(~exist('lmin', 'var'))
                lmin = 2;
            end

            N = length(H);
            
            value = (lmin:N)*H(lmin:N)/sum(H(lmin:N));
        end
        
        function value = LMAX(obj)
            if(isempty(obj.diagonalLinesHistogram))
                notify(obj, 'callHistCalc')
            end
            H = obj.diagonalLinesHistogram;
            
            value = find(H(1:(end-1)), 1, 'last');
        end
        
        function value = DIV(obj)
            value = 1/LMAX(obj);
        end
        
        function value = ENTD(obj, lmin)
            if(isempty(obj.diagonalLinesHistogram))
                notify(obj, 'callHistCalc')
            end
            H = obj.diagonalLinesHistogram;
            
            if(~exist('lmin', 'var'))
                lmin = 2;
            end
            
            H = H(lmin:end);
            
            idx = (H > 0);
            
            p = H(idx)./sum(H);
            
            value = -sum(p.*log(p));
        end
        %}
        
    end %END protected methods
    
end
