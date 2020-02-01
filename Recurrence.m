classdef Recurrence
%RECURRENCE Recurrence
%   This class contains the recurrence parameters that are used to generate the
%   distance matrix and the recurrence plot.
%
% SYNTAX
%   obj = Recurrence(embeddingDimension, timeDelay, normType);
%   obj = Recurrence(embeddingDimension, timeDelay, threshold, normType);
%
% PROPERTIES
%   embeddingDimension
%       Embedding dimension (positive integer greater than or equal to 1)
%
%   timeDelay
%       Time delay (positive integer greater than or equal to 1)
%
%   threshold
%       Threshold (non-negative real number or interval)
%
%   normType
%       Norm used to compute the distance between state-space vectors
%
% METHODS (PUBLIC)
%   Recurrence()
%       Class constructor function
% 
% METHODS (PUBLIC, STATIC)
%   plotr()
%       Plot distance matrix or recurrence plot
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Feb 1, 2020
% ============================================================================ %

    properties
        % Recurrence parameters
        embeddingDimension  % Embedding dimension
        timeDelay           % Time delay
        threshold           % Threshold
        normType            % Norm (L1, L2, L-infinity)
    end % END properties
    
    methods
        %
        % Class constructor function
        %
        
        function obj = Recurrence(varargin)
            switch nargin
                case 3  % Distance matrix
                    obj.embeddingDimension = varargin{1};
                    obj.timeDelay          = varargin{2};
                    obj.normType           = varargin{3};
                case 4  % Recurrence plot
                    obj.embeddingDimension = varargin{1};
                    obj.timeDelay          = varargin{2};
                    obj.threshold          = varargin{3};
                    obj.normType           = varargin{4};
                otherwise
                    ERR_MSG = "Wrong number of input arguments.";
                    error(ERR_MSG);
            end
            
            % Check recurrence parameters
            obj = checkparameters(obj); 
        end % END Recurrence()
        
        %
        % Gets and sets
        % 
        
        function obj = setembeddingdimension(obj, newEmbeddingDimension)
        %SETEMBEDDINGDIMENSION
        % -------------------------------------------------------------------- %
            obj.embeddingDimension = newEmbeddingDimension;
        end % END setembeddingdimension()
        
        function obj = settimedelay(obj, newTimeDelay)
        %SETTIMEDELAY
        % -------------------------------------------------------------------- %
            obj.timeDelay = newTimeDelay;
        end % END settimedelay()
        
        function obj = setthreshold(obj, newThreshold)
        %SETTHRESHOLD
        %
        % -------------------------------------------------------------------- %
            obj = setthreshold@Recurrence(obj, newThreshold);
        end % END setthreshold()
        
        function obj = setnormtype(obj, newNormType)
        %SETNORMTYPE
        % -------------------------------------------------------------------- %
            obj.normType = newNormType;
        end % END setnormtype()
        
    end % END public methods
    
    methods (Access = protected)
        %
        % Checking parameters
        %
        
        function obj = checkparameters(obj)
        %CHECKPARAMETERS Check recurrence parameters
        %   This function validates the embedding dimension, the time delay, the
        %   threshold, and the norm
        % -------------------------------------------------------------------- %
            if(isa(obj, 'DistanceMatrix'))
                obj = checkembeddingdimension(obj);
                obj = checktimedelay(obj);
                obj = checknorm(obj);
            elseif(isa(obj, 'RecurrencePlot'))
                obj = checkembeddingdimension(obj);
                obj = checktimedelay(obj);
                obj = checkthreshold(obj);
                obj = checknorm(obj);
            end
        end % END checkparameters()
        
        function obj = checkembeddingdimension(obj)
        %CHECKEMBEDDINGDIMENSION
        %   This function validates the embedding dimension.
        % -------------------------------------------------------------------- %
            ERR_MSG = "Invalid parameter: embedding dimension.";
        
            if(isempty(obj.embeddingDimension))
                ERR_MSG = strcat(ERR_MSG, " Argument cannot be empty.");
            end
            
            if(ischar(obj.embeddingDimension))
                obj.embeddingDimension = sprintf("%c", obj.embeddingDimension);
            end
        
            if(isstring(obj.embeddingDimension))
                obj.embeddingDimension = str2double(obj.embeddingDimension);
            end
                
            [m,n] = size(obj.embeddingDimension);
            if(m ~= 1 || n ~= 1)
                ERR_MSG = strcat(ERR_MSG, " It must be a positive integer", ...
                    " number greater than or equal to 1.");
                error(ERR_MSG);
            end
            
            if(obj.embeddingDimension < 1)
                ERR_MSG = strcat(ERR_MSG, " It must be a positive integer", ...
                    " number greater than or equal to 1.");
                error(ERR_MSG);
            end
            
            if(mod(obj.embeddingDimension, 1) ~= 0)
                ERR_MSG = strcat(ERR_MSG, " It must be a positive integer", ...
                    " number.");
                error(ERR_MSG);
            end

        end % END checkembeddingdimension()
        
        function obj = checktimedelay(obj)
        %CHECKTIMEDELAY Check time delay
        %   This function validates the time delay, which is integer number
        %   greater than or equal to 0. If the embedding dimension is equal to
        %   1, then the value of the time delay doesn't matter (in this case, we
        %   set time delay to 0).
        %
        %   User can input the time delay value as a number (e.g., 2), as string
        %   (e.g., "2"), or as a char (e.g., '2').
        % -------------------------------------------------------------------- %
            ERR_MSG = "Invalid parameter: time delay.";
            
            if(isempty(obj.timeDelay))
                ERR_MSG = strcat(ERR_MSG, " Argument cannot be empty.");
                error(ERR_MSG);
            end
            
            if(ischar(obj.timeDelay))
                obj.timeDelay = sprintf("%c", obj.timeDelay);
            end
        
            if(isstring(obj.timeDelay))
                obj.timeDelay = str2double(obj.timeDelay);
                if(isempty(obj.timeDelay))
                    error(ERR_MSG);
                end
            end
            
            if(isnumeric(obj.timeDelay))
                % Negative number?
                if(obj.timeDelay < 0)
                    ERR_MSG = strcat(ERR_MSG, " Value must be positive.");
                    error(ERR_MSG);
                end
                
                % Vector?
                [m, n] = size(obj.timeDelay);
                if(m > 1 || n > 1)
                    ERR_MSG = strcat(ERR_MSG, " Value must be a scalar.");
                    error(ERR_MSG);
                end
                
                % Not an integer?
                if(mod(obj.timeDelay, 1) ~= 0)
                    ERR_MSG = strcat(ERR_MSG, " Value must be an integer.");
                    error(ERR_MSG);
                end
            end
            
            if(obj.embeddingDimension == 1 && obj.timeDelay ~= 0)
                obj.timeDelay = 0;
            end
            
            if(obj.embeddingDimension > 1 && obj.timeDelay == 0)
                ERR_MSG = strcat(ERR_MSG, " Value must be greater than", ...
                    " or equal to 1.");
                error(ERR_MSG);
            end
        end % END checktimedelay()
        
        function obj = checkthreshold(obj)
        %CHECKTHRESHOLD Check the threshold parameter
        %   
        % -------------------------------------------------------------------- %
            ERR_MSG = "Invalid parameter: threshold.";
            
            if(isempty(obj.threshold))
                ERR_MSG = strcat(ERR_MSG, " Argument cannot be empty.");
                error(ERR_MSG);
            end
        
            if(ischar(obj.threshold))
                obj.threshold = sprintf("%c", obj.threshold);
            end
            
            if(isstring(obj.threshold))
                obj.threshold = str2double(obj.threshold);
            end
            
            [m, n] = size(obj.threshold);
            if(m == 1 && n == 1)
                obj.threshold(2) = 0;
                [m, n] = size(obj.threshold);
            end
            
            if(m < 0 || m > 2 || n < 0 || n > 2)
                ERR_MSG = strcat(ERR_MSG, " It must be a single real", ...
                    " number or an interval (e.g., [0.1, 0.25]).");
                error(ERR_MSG);
            end
            
            if(m == 2 && n == 2)
                ERR_MSG = strcat(ERR_MSG, " It must be a single real", ...
                    " number or an interval (e.g., [0.1, 0.25]).");
                error(ERR_MSG);
            end
            
            if((m == 2 && n == 1) || (m == 1 && n == 2))
                % Corridor thresholded version of the recurrence plot
                obj.threshold = sort(obj.threshold, 'ascend');
            end            
        end % END validatethreshold()
        
        function obj = checknorm(obj)
        %CHECKNORM Check the norm used to calculate the distance between
        %state-space vectors
        %   This function validates the norm. It must be choosen among L1
        %   (a.k.a. taxicab norm), L2 (a.k.a. Euclidean norm), and L-infinity
        %   (a.k.a. maximum norm).
        %   
        %   User can use these norms as follows:
        %       * as strings: "L1", "L2", and "L3";
        %       * as chars: 'L1', 'L2', and 'L3';
        %       * as numbers: 1, 2, and 3 (respectively to L1, L2, and L3).
        % -------------------------------------------------------------------- %
            ERR_MSG = "Invalid parameter: norm.";
            
            if(isstring(obj.normType))
                obj.normType = char(obj.normType);
            end
            
            if(isnumeric(obj.normType))
                if(mod(obj.normType, 1) == 0)
                    if(obj.normType < 1 || obj.normType > 3)
                        ERR_MSG = strcat(ERR_MSG, " Norm can be an integer ", ...
                            "number between 1 and 3: 1 = L1, 2 = L2, and ", ...
                            "3 = L-infinity.");
                        error(ERR_MSG);
                    else
                        switch obj.normType
                            case 1
                                obj.normType = 'L1';
                            case 2
                                obj.normType = 'L2';
                            case 3
                                obj.normType = 'L-infinity';
                        end
                    end
                else
                    ERR_MSG = strcat(ERR_MSG, " Norm can be an integer ", ...
                    "number: 1 = L1,2 = L2, 3 = L-infinity.");
                    error(ERR_MSG);
                end
                
            end
            
            VALID_NORMS = {'L1','l1', 'taxicab', ...
                'L2', 'l2', 'euclidean', 'Euclidean', ...
                'L-infinity', 'l-infinity', 'Linf', 'linf', 'Chebychev'};
            if(~ismember(obj.normType, VALID_NORMS))
                ERR_MSG = strcat(ERR_MSG, " Choose among L1, L2, and ", ...
                    "L-infinity.");
                error(ERR_MSG);
            end

            % Standardize
            switch obj.normType
                case {'l1', 'taxicab'}
                    obj.normType = 'L1';
                case {'l2', 'euclidean', 'Euclidean'}
                    obj.normType = 'L2';
                case {'l-infinity', 'Linf', 'linf', 'Chebychev'}
                    obj.normType = 'L-infinity';
            end
        end % END checknorm()
                
    end % END protected methods

    methods (Static)
        %
        % Plot
        %
        
        function plotr(M)
        %PLOTR Plot recurrence object
        %   Plots the self-distance / cross-distance matrix or the recurrence /
        %   cross-recurrence plot. If the object is a recurrence /
        %   cross-recurrence plot, the image is set to black and white.
        % -------------------------------------------------------------------- %
            imagesc(M)
            
            if(islogical(M))
                % Set colors to B&W
                colormap([1, 1, 1; 0.5*rand(1, 3)])
            else
                colormap('parula')
                colorbar
            end

            [m, n] = size(M);
            if(m == n)  % symmetric matrix
                axis square
            end

            % Axis
            [m, n] = size(M);
            ax = [1, round(m/2), m];
            ay = [1, round(n/2), n];
            set(gca, 'XTick', ax, 'XTickLabel', ax)
            set(gca, 'YTick', ay, 'YTickLabel', ay)
        end % END plotr()
        
    end
    
end
