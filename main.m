%% Recurrence Analysis
%
% File: main.m
%
% The purpose of this file is to demonstrate how to use the following classes:
%   RP,
%   RQA, and
%   RNA.
%
% CONTACT
%   Patrick Franco Coutinho
%   pfcoutinho@outlook.com
%
% Last update: Feb 17, 2020
%% Example #1
% Stochastic process (uniform distribution).

close all
clear
clc

%
% Time series data
%
data = rand(100, 1);

%
% Recurrence parameters
%
embeddingDimension = 1;
timeDelay          = 1;
threshold          = 0.1;
normType           = 'L1';

%
% Recurrence plot
%
try
    rp = RP(embeddingDimension, timeDelay, threshold, normType, data);
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
xlim([1, numel(data)])
set(gca, "XTick", [1, round(numel(data)/(2)), numel(data)]);

subplot(5, 1, [2, 5])
rp.plotr()

%% Example #2
% Periodic time series (sinusoidal signal).

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
    rp = RP(embeddingDimension, timeDelay, threshold, normType, data);
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
rp.plotr()

%% Example #3
% Auto-regressive model.

close all
clear
clc

%
% Time series data
%
data    = zeros(100, 1);
data(1) = 0.1;
for i = 2:100
    data(i) = data(i-1) + 0.3*randn();
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
    rp = RP(embeddingDimension, timeDelay, threshold, normType, data);
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
rp.plotr()

%% Example #4
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
threshold          = 0.4;
normType           = 'L-infinity';

%
% Recurrence plot
%
try
    rp = RP(embeddingDimension, timeDelay, threshold, normType, data);
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
rp.plotr()

%% Example #5
% Stochastic process and corridor threshold.

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
threshold          = [0.1, 0.3];
normType           = 'L1';

%
% Recurrence plot
%
try
    rp = RP(embeddingDimension, timeDelay, threshold, normType, data);
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
rp.plotr()

%% Example #6
% Sthocastic processes (uniform distribution) and cross recurrence plot.

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
threshold          = 0.1;
normType           = 'L1';

%
% Recurrence plot
%
try
    rp = RP(embeddingDimension, timeDelay, threshold, normType, data1, data2);
catch ERR
    error(ERR.message)
end

%
% Plotting the recurrence matrix
%
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
rp.plotr()

set(gcf, 'Position', [400, 140, 560, 530]);

%%
