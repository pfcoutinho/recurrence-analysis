%% Recurrence Analysis
%
% File: main.m
%
% The purpose of this file is to demonstrate how to use the following classes:
%   DistanceMatrix,
%   RecurrencePlot,
%   RecurrenceQuantificationAnalysis, and
%   RecurrenceNetworkAnalysis.
%
% AUTHOR
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: May 11, 2020
%% Example #1
% Distance matrix of a stochastic process (drawn from the uniform distribution)

close all
clear
clc

%
% Time series data
%
data = rand(30, 1);

%
% Parameters
%
embeddingDimension = 1;
timeDelay          = 1;
normType           = 'l1';

%
% Distance matrix
%
try
    dm = DistanceMatrix(embeddingDimension, timeDelay, normType, data);
catch ERR
    error(ERR.message)
end

% Plot the distance matrix
figure(1)

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, numel(data)])
set(gca, "XTick", [1, round(numel(data)/(2)), numel(data)]);

subplot(5, 1, [2, 5])
dm.plot()

% Setting a new embedding dimension
dm.embeddingDimension = 3;

% Plot the distance matrix
figure(2)
subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, numel(data)])
set(gca, "XTick", [1, round(numel(data)/(2)), numel(data)]);

subplot(5, 1, [2, 5])
dm.plot()

% Setting a new norm
dm.normType = 'l2';

% Plot the distance matrix
figure(3)
subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, numel(data)])
set(gca, "XTick", [1, round(numel(data)/(2)), numel(data)]);

subplot(5, 1, [2, 5])
dm.plot()

%% Example #2
% Recurrence plot of a stochastic process (drawn from the uniform distribution)

close all
clear
clc

%
% Time series data
%
data = rand(50, 1);

%
% Recurrence parameters
%
embeddingDimension = 1;
timeDelay          = 1;
threshold          = 0.646;
normType           = 'L1';

%
% Recurrence plot
%
try
    RPobj = RecurrencePlot(embeddingDimension, timeDelay, ...
        'recurrence rate', threshold, normType, data)
catch ERR
    error(ERR.message)
end

% Plot the recurrence matrix
figure(1)

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, numel(data)])
set(gca, "XTick", [1, round(numel(data)/(2)), numel(data)]);

subplot(5, 1, [2, 5])
RPobj.plot('color')

%
% setting a new threshold
%
RPobj.threshold = 0.08

% Plot the recurrence matrix
figure(2)

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, numel(data)])
set(gca, "XTick", [1, round(numel(data)/(2)), numel(data)]);

subplot(5, 1, [2, 5])
RPobj.plot('color')

%
% setting a new recurrence rate
%
RPobj.recurrenceRate = RPobj.recurrenceRate+eps;

% Plot the recurrence matrix
figure(3)

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, numel(data)])
set(gca, "XTick", [1, round(numel(data)/(2)), numel(data)]);

subplot(5, 1, [2, 5])
RPobj.plot('color')

%
% setting a new time series
%
newData = rand(40, 1);
RPobj.data = newData;

% Plot the recurrence matrix
figure(4)

subplot(5, 1, 1)
plot(newData)
ylabel("x(t)")
xlabel("t")
xlim([1, numel(newData)])
set(gca, "XTick", [1, round(numel(newData)/(2)), numel(newData)]);

subplot(5, 1, [2, 5])
RPobj.plot()

%
% setting a new embedding dimension
%
RPobj.embeddingDimension = 2;

% Plot the recurrence matrix
figure(5)

subplot(5, 1, 1)
plot(newData)
ylabel("x(t)")
xlabel("t")
xlim([1, numel(newData)])
set(gca, "XTick", [1, round(numel(newData)/(2)), numel(newData)]);

subplot(5, 1, [2, 5])
RPobj.plot()

%% Example #3
% Periodic time series (sinusoidal signal)

close all
clear
clc

%
% Time series data
%
data = sin(2*pi*0.05*(1:100));

%
% Recurrence parameters
%
embeddingDimension = 2;
timeDelay          = 3;
threshold          = 0.2;
normType           = 'L2';

%
% Recurrence plot
%
try
    RPobj = RecurrencePlot(embeddingDimension, timeDelay, ...
                'threhsold', threshold, normType, data);
catch ERR
    error(ERR.message)
end

%
% Plotting the recurrence matrix
%
figure(1)

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RPobj.plot()

%% Example #4
% Auto-regressive model

close all
clear
clc

%
% Time series data
%
data    = zeros(100, 1);
data(1) = 0.1;
for i = 2:100
    data(i) = data(i-1) + 0.2*randn();
end

%
% Recurrence parameters
%
embeddingDimension = 2;
timeDelay          = 1;
threshold          = 0.3;
normType           = 'L2';

%
% Recurrence plot
%
try
    RPobj = RecurrencePlot(embeddingDimension, timeDelay, ...
                'threshold', threshold, normType, data);
catch ERR
    error(ERR.message)
end

%
% Plotting the recurrence matrix
%
figure()

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RPobj.plot('color')

%% Example #5
% Chaotic time series (logistic map)

close all
clear
clc

%
% Time series data
%
data    = zeros(100, 1);
data(1) = 0.2;
for i = 2:100
    data(i) = 4*data(i-1)*(1-data(i-1));
end

%
% Recurrence parameters
%
embeddingDimension = 3;
timeDelay          = 1;
threshold          = 0.2;
normType           = 'L-infinity';

%
% Recurrence plot
%
try
    RPobj = RecurrencePlot(embeddingDimension, timeDelay, ...
                'threshold', threshold, normType, data);
catch ERR
    error(ERR.message)
end

%
% Plotting the recurrence matrix
%
figure(1)

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RPobj.plot()

%% Example #6
% Stochastic process (drawn from the Gaussian distribution) and corridor threshold

close all
clear
clc

%
% Time series data
%
data = randn(100, 1);

%
% Recurrence parameters
%
embeddingDimension = 1;
timeDelay          = 1;
threshold          = [0.1, 0.3];    % corridor threshold
normType           = 'L1';

%
% Recurrence plot
%
try
    RPobj = RecurrencePlot(embeddingDimension, timeDelay, ...
                'recurrence rate', threshold, normType, data);
catch ERR
    error(ERR.message)
end

%
% Plotting the recurrence matrix
%
figure(1)

subplot(5, 1, 1)
plot(data)
ylabel("x(t)")
xlabel("t")
xlim([1, 100])
set(gca, "XTick", [1, 10:10:100])

subplot(5, 1, [2, 5])
RPobj.plot()

%% Example #7
% Sthocastic processes (drawn from the uniform distribution) and cross 
% recurrence plot

close all
clear
clc

%
% Time series
%
data1 = rand(100, 1);
data2 = rand(150, 1);

%
% Recurrence parameters
%
embeddingDimension = 1;
timeDelay          = 1;
threshold          = 0.5;
normType           = 'L1';

%
% Recurrence plot
%
try
    RPobj = RecurrencePlot(embeddingDimension, timeDelay, ...
                'threshold', threshold, normType, data1, data2)
catch ERR
    error(ERR.message)
end

%
% Plotting the recurrence matrix
%
figure(1)

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
RPobj.plot('color')

set(gcf, 'Position', [300, 80, 530, 600]);
