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

x = sin(2*pi*0.05*(1:100));

RP = recplot(2, 3, 0.2, 'L1', x);

figure()

subplot(5, 1, 1)
plot(x)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
plotrp(RP)

%% Distance matrix and recurrence plot for a random time series
close all
clear
clc

x = rand(100, 1);

RP = recplot(1, 1, 0.2, 'L1', x);

figure()

subplot(5, 1, 1)
plot(x)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
plotrp(RP)

%% Corridor thresholded version of the recurrence plot for one time series
close all
clear
clc

x = rand(100, 1);

RP = recplot(1, 1, [0.1, 0.2], 'Linf', x);

figure()

subplot(5, 1, 1)
plot(x)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
plotrp(RP)

%% Cross-recurrence plot between two time series
close all
clear
clc

x = rand(100, 1);
y = rand(150, 1);

XRP = recplot(1, 1, 0.2, 'Linf', x, y);

figure()
plotrp(XRP)

%% Corridor thresholded version of cross-recurrence plot between two time series
close all
clear
clc

x = rand(100, 1);
y = rand(150, 1);

XRP = recplot(1, 1, [0.1, 0.2], 'Linf', x, y);

figure()
plotrp(XRP)
