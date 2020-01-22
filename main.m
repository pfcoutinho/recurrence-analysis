%
% The purpose of this file is to test the following functions:
%   distmatrix.m
%   recplot.m
%   plotrp.m
%

%% Distance matrix and recurrence plot for one time series
close all
clear
clc

X = rand(100, 1);

DMofX = distmatrix(1, 1, 'Linf', X);
RPofX = recplot(1, 1, 0.2, 'Linf', X);

figure()
plotrp(DMofX)


figure()
plotrp(RPofX)

%% Corridor thresholded version of the recurrence plot for one time series
close all
clear
clc

X = rand(100, 1);

RPofX = recplot(1, 1, [0.0, 0.2], 'Linf', X);

figure()
plotrp(RPofX, 'bw')
plotrp(RPofX)

%% Cross-recurrence plot between two time series
close all
clear
clc

X = rand(100, 1);
Y = rand(150, 1);

XRPofXandY = recplot(1, 1, 0.2, 'Linf', X, Y);

figure()
plotrp(XRPofXandY)
