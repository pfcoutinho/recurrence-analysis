% main.m
%
% The purpose of this file is to test the following classes:
%   DistanceMatrix
%   RecurrencePlot
%   RecurrenceAnalysis
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 28, 2020
%% Distance matrix and recurrence plot for a sinusoidal time series
close all
clear
clc

data = sin(2*pi*0.05*(1:100));

RP = RecurrencePlot(2, 3, 0.2, 'L2', data);
D  = DistanceMatrix(2, 3, 'L2', data);

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RP.plot()

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
D.plot()

%% Distance matrix and recurrence plot for a random time series
close all
clear
clc

data = rand(100, 1);

RP = RecurrencePlot(1, 1, 0.2, 'L1', data);

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RP.plot()

%% Corridor thresholded version of the recurrence plot for one time series
close all
clear
clc

data = rand(100, 1);

RP = RecurrencePlot(1, 1, [0.1, 0.2], 'L-infinity', data);

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RP.plot()

%% Cross-recurrence plot between two time series
%{
close all
clear
clc

data1 = rand(100, 1);
data2 = rand(150, 1);

XRP = (1, 1, 0.2, 'L-infinity', data1, data2);

figure()
plotrp(XRP)
%}

%% Corridor thresholded version of cross-recurrence plot between two time series
%{
close all
clear
clc

data1 = rand(100, 1);
data2 = randn(150, 1);

XRP = recplot(1, 1, [0.1, 0.2], 'Linf', data1, data2);

figure()
plotrp(XRP)
%}
