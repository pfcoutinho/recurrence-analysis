classdef RecurrenceQuantificationAnalysis < handle
%RQA Recurrence Quantification Analysis
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
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Mar 10, 2020
% ============================================================================ %
    %
    % Properties
    %
    properties
        % Histograms
        diagonalLinesHistogram  % Histogram of diagonal lines
        verticalLinesHistogram  % Histogram of vertical lines
    end
    
    properties (Abstract = true)
        M                       % Recurrence matrix
        timeSeries              % Time series
    end
    
    %
    % Events
    %
    events
        callHistCalc
    end
    
    %
    % Methods
    %
    methods
        function self = RecurrenceQuantificationAnalysis()
            % 
            addlistener(self, 'callHistCalc', @self.histevthandle);
        end
        
        %
        % Histograms
        %
        function value = DHIST(self)
        %DHIST Histogram of diagonal lines
        % -------------------------------------------------------------------- %
            if(isempty(self.diagonalLinesHistogram))
                notify(self, 'callHistCalc');
            end
            value = self.diagonalLinesHistogram;
        end % END DHIST()
        
        
        function value = VHIST(self)
        %VHIST Histogram of vertical lines
        % -------------------------------------------------------------------- %
            value = self.verticalLinesHistogram;
        end % END VHIST()
        
        %
        % Measures based on the density of recurrence points
        %
        function value = RR(obj)
        %RR Recurrence Rate
        % -------------------------------------------------------------------- %
            if(isempty(obj.M))
                error("Cannot calculate the recurrence rate (RR).")
            end
        
            [m, n] = size(obj.M);
            value = nnz(obj.M)/(m*n);
        end % END RR()
        
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
        
        %
        % Adaptive threshold
        %
        %{
        function [threshold] = adptthreshold(data, rate)
            
        end
        %}
    end
    
    methods (Static = true)

    end
    
    methods (Access = private)
        function self = histevthandle(self, ~, ~)
            self = chkhist(self, 'diagonal');
            self = chkhist(self, 'vertical');
        end
        
        function self = chkhist(self, opt)
        % -------------------------------------------------------------------- %
            switch opt
                case 'diagonal'
                    if(~hashist(self, opt))
                        self.diagonalLinesHistogram = hod(self);
                    end
                case 'vertical'
                    if(~hashist(self, opt))
                        self.verticalLinesHistogram = hov(self);
                    end
            end
            
            function value = hashist(self, opt)
                switch opt
                    case 'diagonal'
                        if(isempty(self.diagonalLinesHistogram))
                            value = false;
                        else
                            value = true;
                        end
                    case 'vertical'
                        if(isempty(self.verticalLinesHistogram))
                            value = false;
                        else
                            value = true;
                        end
                end
            end
        end % END chkhist()
        
        function H = hod(obj)
        %HOD Histogram of diagonal lines
        % -------------------------------------------------------------------- %
            [m, n] = size(obj.M);
            
            H = zeros(max(m, n), 1);
            
            for i = -(n-1):1:(n-1)
                diagonalLine = [0; find(~diag(obj.M, i)); numel(diag(obj.M, i)) + 1];
                lineLengths  = diff(diagonalLine) - 1;
                idx          = (lineLengths > 0);
                H(lineLengths(idx)) = H(lineLengths(idx)) + 1;
            end
        end % END hod()
        
        function H = hov(obj)
        %HOD Histogram of vertical lines
        % -------------------------------------------------------------------- %
            [m, n] = size(obj.M);
            
            H = zeros(max(m, n), 1);
            
            for i = 1:m
                verticalLine = [0; find(~obj.M(1:end, i)); numel(obj.M(1:end, i)) + 1];
                lineLengths  = diff(verticalLine) - 1;
                idx          = (lineLengths > 0);
                H(lineLengths(idx)) = H(lineLengths(idx)) + 1;
            end
        end % END hov()
    end
    
end
