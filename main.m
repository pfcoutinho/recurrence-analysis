%% Recurrence Analysis
%
% File: main.m
%
% The purpose of this file is to demonstrate how to use the following classes:
%   DistanceMatrix,
%   RecurrencePlot, and
%   RecurrenceAnalysis.
%
% *CONTACT*
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Jan 28, 2020
%% HOW TO USE
%
% To run one section at a time, press CTRL + ENTER.
%
% ## Distance Matrix
%
% To generate a distance matrix from a single time series,
%
%   D = DistanceMatrix(embeddingDimension, timeDelay, normType, data)
%
% where:
%
%   * embeddingDimension is an integer number greater than or equal to 1;
%   * timeDelay is an integer number greater than or equal to 0;
%   * normType is one of the norms L1 (a.k.a., taxicab norm), L2 (a.k.a., 
%       Euclidean norm), and L-infinity (a.k.a., maximum norm);
%   * data is the time series.
%
% The distance matrix can be accessed by the command D.M;
%
% To calculate the distance matrix between two time series:
%
%   D = DistanceMatrix(embeddingDimension, timeDelay, normType, data_x, data_y)
%
% ## Recurrence plot
%
%
%
%
%% Distance matrix and recurrence plot for a sinusoidal time series
close all
clear
clc

%
% Time series data
%
data = sin(2*pi*0.05*(1:100)');

%
% Distance matrix
%
try
    D = DistanceMatrix(2, 3, 'L2', data);
catch ERR
    error(ERR.message)
end

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
D.plotr()   % or plot(D);

%
% Recurrence plot
%
try
    RP = RecurrencePlot(2, 3, 0.2, 'l2', data);
catch ERR
    error(ERR.message)
end

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RP.plotr()   % or plotr(RP);

%% Distance matrix and recurrence plot for a random time series
close all
clear
clc

%
% Time series
%
data = rand(100, 1);

%
% Distance matrix
%
try
    D = DistanceMatrix(1, 1, 'l1', data);
catch

end

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
D.plotr()

%
% Recurrence plot
%
try
    RP = RecurrencePlot(1, 1, 0.2, 'l1', data);
catch

end

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RP.plotr()

%% Corridor thresholded version of the recurrence plot
close all
clear
clc

%
% Time series
%
data = rand(100, 1);

%
% Recurrence plot
%
try
    RP = RecurrencePlot(1, 1, [0.1, 0.2], 'L-infinity', data);
catch ERR
    error(ERR.message)
end

figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RP.plotr()

%% Cross-distance matrix and cross-recurrence plot between two time series

close all
clear
clc

%
% Time series
%
data1 = rand(100, 1);
data2 = rand(150, 1);

%
% Distance matrix
%
try
    D = DistanceMatrix(2, 3, 'L2', data1, data1);
catch ERR
    error(ERR.message)
end

figure()

subplot(6, 1, 1)
plot(data1)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(6, 1, 2)
plot(data2)
ylabel("y(t)")
xlabel("t")
xlim([1, 150])
set(gca, "XTick", [1, 10:10:150])

subplot(6, 1, [3, 6])
D.plotr()

set(gcf, 'Position', [400, 140, 560, 530]);

%
% Recurrence plot
%
try
    RP = RecurrencePlot(2, 3, 0.2, 'L2', data1, data1);
catch ERR
    error(ERR.message)
end

figure()

subplot(6, 1, 1)
plot(data1)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(6, 1, 2)
plot(data2)
ylabel("y(t)")
xlabel("t")
xlim([1, 150])
set(gca, "XTick", [1, 10:10:150])

subplot(6, 1, [3, 6])
RP.plotr()

set(gcf, 'Position', [400, 140, 560, 530]);

%% Cross-distance matrix and cross-recurrence plot between two time series (corridor version)

close all
clear
clc

%
% Time series
%
data1 = rand(100, 1);
data2 = rand(150, 1);

%
% Distance matrix
%
try
    D = DistanceMatrix(2, 3, 'L2', data1, data1);
catch ERR
    error(ERR.message)
end

figure()

subplot(6, 1, 1)
plot(data1)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(6, 1, 2)
plot(data2)
ylabel("y(t)")
xlabel("t")
xlim([1, 150])
set(gca, "XTick", [1, 10:10:150])

subplot(6, 1, [3, 6])
D.plotr()

set(gcf, 'Position', [400, 140, 560, 530]);

%
% Recurrence plot
%
try
    RP = RecurrencePlot(2, 3, [0.2, 0.3], 'L2', data1, data1);
catch ERR
    error(ERR.message)
end

figure()

subplot(6, 1, 1)
plot(data1)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(6, 1, 2)
plot(data2)
ylabel("y(t)")
xlabel("t")
xlim([1, 150])
set(gca, "XTick", [1, 10:10:150])

subplot(6, 1, [3, 6])
RP.plotr()

set(gcf, 'Position', [400, 140, 560, 530]);
