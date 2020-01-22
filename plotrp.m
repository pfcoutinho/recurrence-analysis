function plotrp(M)
%PLOTRP Plot the recurrence plot or the cross recurrence plot
%   This function plots the recurrence plot of a data series or the cross
%   recurrence plot of two data series.
%
% SYNTAX
%   plotrp(M)
%
% INPUT:
%	M - a matrix
%
% DEPENDENCIES:
%   None.
% 
% REFERENCES:
%   None.
%
% CONTACT:
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%   Last update: Jan 22, 2020
% ============================================================================ %

%     fh = gcf;
    
    if islogical(M)
        imagesc(M)
        colormap([1, 1, 1; 0, 0, 0]) 
    else
        imagesc(M)
        colormap pink
        colorbar
    end
    
    if issymmetric(double(M))
        axis square
    end

    % Removes x and y labels
    set(gca,'XTickLabel',[],'YTickLabel',[])
    
end
