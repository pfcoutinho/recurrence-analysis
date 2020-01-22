function [RP] = recplot(varargin)
%RECPLOT Recurrence plot
%   The function computes the recurrence plot of a time series or the cross
%   recurrence plot between two time series.
%
% SYNTAX:
%   RP  = recplot(x, de, tau, threshold, dfunc)
%   XRP = recplot(x, y, de, tau, threshold, dfunc)
%
% INPUT:
%   EMBEDDING_DIMENSION - embedding dimension
%   TIME_DELAY          - time delay
%   THRESHOLD           - threshold (scalar, or 1 x 2 or 2 x 1 vector)
%   DTYPE               - similarity or distance function
%   x, y                - time series (vectors)
%
% OUTPUT:
%   RP        - recurrence or cross-recurrence plot (binary M x M matrix)
%
% DEPENDENCIES:
%   distancematrix.m
%
% REFERENCES:
%   Site: <http://www.agnld.uni-potsdam.de/~marwan/matlab-tutorials/html/rp.html>
%   N. Marwan, M. C. Romano, M. Thiel, and J. Kurths, Recurrence Plots for
%       the Analysis of Complex Systems, Physics Reports, Vol. 438, pp.237-329,
%       2007 
%
% CONTACT:
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%   Last update: Jan 22, 2020
% ============================================================================ %

    %
    % Checking the number of input and output arguments
    %
    
    % Input
    narginchk(5, 6);
    
    % Output
    nargoutchk(1, 1);
    % ------------------------------------------------------------------------ %
    
    %
    % Getting the parameters of the recurrence plot and the time series
    %
    
    % Recurrence plot parameters
    EMBEDDING_DIMENSION = varargin{1};
    TIME_DELAY = varargin{2};
    THRESHOLD  = varargin{3};
    DTYPE      = varargin{4};
    
    % Time series
    x = varargin{5};
    if nargin == 6
        y = varargin{6};
    end
    % ------------------------------------------------------------------------ %
    
    %
    % Checking the arguments
    %
    
    % Threshold
    if ~isnumeric(THRESHOLD)
        error('Threshold must be numeric');
    end
    
    if ~isempty(THRESHOLD)
        if isnumeric(THRESHOLD)
            if isscalar(THRESHOLD)
                FLAG_CORRIDOR = false;
            else
                % corridor threhsold;
                if (size(THRESHOLD, 1) == 1 && size(THRESHOLD,2) == 2) || ...
                        (size(THRESHOLD, 1) == 2 && size(THRESHOLD,2) == 1)
                    if size(THRESHOLD,1) == 2
                        THRESHOLD = THRESHOLD';
                    end                        
                    FLAG_CORRIDOR = true;
                else
                    error('Threshold must be 1 x 2 or 2 x 1 vector.');
                end
            end
        else
            error('Threshold (epsilon) must be numeric');
        end
    end
    
% ============================================================================ %

    % Obtain the distance matrix
    switch nargin
        case 5  % recurrence plot;
            DM = distmatrix(EMBEDDING_DIMENSION, TIME_DELAY, DTYPE, x);
        case 6  % cross recurrence plot;
            DM = distmatrix(EMBEDDING_DIMENSION, TIME_DELAY, DTYPE, x, y);
    end
    
    % Apply the trheshold
    if ~FLAG_CORRIDOR
        % thresholded version;
        RP = (abs(DM) <= THRESHOLD);
    else
        % corridor thresholded version;
        RP_min = (abs(DM) >= THRESHOLD(1));
        RP_max = (abs(DM) <= THRESHOLD(2));
        RP = and(RP_min, RP_max);
    end
        
end

