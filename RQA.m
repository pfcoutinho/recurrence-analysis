classdef RQA
%RQA Recurrence Quantification Analysis
%   This class performs recurrence quantification analysis (RQA) of recurrence 
%   plots, which are binary matrices (i.e., A(i, j) = {0, 1}).
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Feb 15, 2020
% ============================================================================ %
    properties (Abstract = true)
        M                   % Recurrence matrix
    end

    properties (Dependent)
        % Based on the density of recurrence points
        RR                  % Recurrence rate
        
        % Histograms
        HDL                 % Histogram of diagonal lines
        HVL                 % Histogram of vertical lines

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
        function value = get.HDL(obj)
            value = RQA.hod(obj.M);
        end
        
        function value = get.HVL(obj)
            value = RQA.hov(obj.M);
        end
        
        function value = get.RR(obj)
            [m, n] = size(obj.M);
            % Discard the LOI?
            if(m == n && sum(diag(obj.M)) == m)
                value = (nnz(obj.M) - m)/(m*n);
            else
                value = nnz(obj.M)/(m*n);
            end
        end
    end
    
    methods (Static = true, Access = private)
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
        end % END hod()
        
        function H = hov(RP)
        %HOD Histogram of vertical lines
        % -------------------------------------------------------------------- %
            [m, n] = size(RP);
            
            H = zeros(max(m, n), 1);
            
            for i = 1:m
                verticalLine = [0; find(~RP(1:end, i)); numel(RP(1:end, i)) + 1];
                lineLengths  = diff(verticalLine) - 1;
                idx          = (lineLengths > 0);
                H(lineLengths(idx)) = H(lineLengths(idx)) + 1;
            end
        end % END hov()
        
        %
        % 
        %
        
%         function value = calcdet(RP, minimumLineLength)
%         end
        
%         function value = calcentd(RP, minimumLineLength)
%         end
        
        %
        %
        %
        
%         function value = calcentv(RP, minimmLineLength)
%         end
        
        %
        % Adaptive threshold
        %
        
%         function adaptivethreshold(data, rate)
%             
%         end
    end
end
