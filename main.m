%
% The purpose of this file is to test the following functions:
%   distmatrix.m
%   recplot.m
%   plotrp.m
%

%% Distance matrix and recurrence plot for a sinusoidal time series
close all
clear
clc

X = sin(2*pi*0.05*(1:100));

DMofX = distmatrix(2, 3, 'L1', X);
RPofX = recplot(2, 3, 0.2, 'L1', X);

figure()
plot(X)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])
set(gcf, "Position", [400, 250, 600, 120]);

figure()
plotrp(DMofX)

figure()
plotrp(RPofX)

%% Distance matrix and recurrence plot for a random time series
close all
clear
clc

X = rand(100, 1);

DMofX = distmatrix(1, 1, 'L1', X);
RPofX = recplot(1, 1, 0.2, 'L1', X);

figure()
plot(X)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])
set(gcf, "Position", [400, 250, 600, 120]);

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
