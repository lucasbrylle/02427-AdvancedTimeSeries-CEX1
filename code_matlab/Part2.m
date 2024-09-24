clc; clear; close all;

% Parameters
T = 1000; % Number of observations
threshold = 0.5; % Threshold for the SETAR model
x_setar = zeros(T, 1); % Initialize the time series
epsilon = randn(T, 1); % Noise term

% SETAR(2,1,1) Model Simulation
for t = 2:T
    if x_setar(t-1) <= threshold
        x_setar(t) = 4 + 0.7 * x_setar(t-1) + epsilon(t);                       % Regime 1
    else
        x_setar(t) = -4 + 0.7 * x_setar(t-1) + epsilon(t);                      % Regime 2
    end
end

% Theoretical Conditional Mean M(x)
M_theoretical = zeros(T-1, 1);
for t = 1:(T-1)
    if x_setar(t) <= threshold
        M_theoretical(t) = 4 + 0.7 * x_setar(t);    % Conditional mean in Regime 1
    else
        M_theoretical(t) = -4 + 0.7 * x_setar(t);   % Conditional mean in Regime 2
    end
end

% Prepare data for Nadaraya-Watson estimation
x_t = x_setar(1:end-1);
y_t_plus_1 = x_setar(2:end);

% Bandwidths to explore
bandwidths = [0.1, 0.5, 1];

% Initialize containers for the estimated M(x) values
M_estimated_gaussian = zeros(T-1, length(bandwidths));
M_estimated_epanechnikov = zeros(T-1, length(bandwidths));

% Gaussian Kernel Function
gaussian_kernel = @(u) exp(-0.5 * u.^2) / sqrt(2 * pi);

% Epanechnikov Kernel Function
epanechnikov_kernel = @(u) 0.75 * (1 - u.^2).* (abs(u) <= 1);

% Nadaraya-Watson Estimation using Gaussian Kernel
for j = 1:length(bandwidths)
    h = bandwidths(j);
    for i = 1:length(x_t)
        u_gaussian = (x_t - x_t(i)) / h;
        weights_gaussian = gaussian_kernel(u_gaussian);
        M_estimated_gaussian(i, j) = sum(weights_gaussian .* y_t_plus_1) / sum(weights_gaussian);
    end
end

% Nadaraya-Watson Estimation using Epanechnikov Kernel
for j = 1:length(bandwidths)
    h = bandwidths(j);
    for i = 1:length(x_t)
        u_epanechnikov = (x_t - x_t(i)) / h;
        weights_epanechnikov = epanechnikov_kernel(u_epanechnikov);
        M_estimated_epanechnikov(i, j) = sum(weights_epanechnikov .* y_t_plus_1) / sum(weights_epanechnikov);
    end
end


% Plotting

% Define the custom color
customRed = [195/255, 31/255, 51/255];
customBlue = [24/255, 54/255, 104/255];

% Font size and weight
fontSize = 13;
fontWeight = 'bold';


% Time vector
time = (1:T-1)';

% 1. Scatter plot of X_t against X_{t-1}
figure
subplot(2,1,1)
scatter(time, x_t, 20, customBlue, 'filled')
xlabel('t', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('X_t', 'FontWeight', fontWeight, 'FontSize', fontSize);
title('SETAR(2,1,1) model');
subplot(2,1,2)
scatter(x_t, y_t_plus_1, 20, customBlue, 'filled');
xlabel('X_{t-1}', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('X_t', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;


% 2. Plot estimated conditional mean with Gaussian kernel
figure
scatter(x_t, y_t_plus_1, 20, customBlue);
hold on;
scatter(x_t, M_theoretical, 20, 'filled');
for j = 1:length(bandwidths)
    scatter(x_t, M_estimated_gaussian(:, j), 20, 'filled');
end
title('Estimated M(X) with Gaussian Kernel');
xlabel('X_t', 'FontWeight',fontWeight,'FontSize',fontSize);
ylabel('M(X)', 'FontWeight',fontWeight,'FontSize',fontSize);
legend(['SETAR simulation', 'Theoretical M(x)', arrayfun(@(h) ['h=', num2str(h)], bandwidths, 'UniformOutput', false)],'Location','southwest');
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

% 3. Plot estimated conditional mean with Epanechnikov kernel
figure
scatter(x_t, y_t_plus_1, 20, customBlue)
hold on;
scatter(x_t, M_theoretical, 20, 'filled');
for j = 1:length(bandwidths)
    scatter(x_t, M_estimated_epanechnikov(:, j), 20, 'filled');
end
title('Estimated M(X) with Epanechnikov Kernel');
xlabel('X_t', 'FontWeight',fontWeight,'FontSize',fontSize);
ylabel('M(X)', 'FontWeight',fontWeight,'FontSize',fontSize);
legend(['SETAR simulation', 'Theoretical M(x)', arrayfun(@(h) ['h=', num2str(h)], bandwidths, 'UniformOutput', false)],'Location','southwest');
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

current_file = 'Part2.m';
new_directory = '/home/olivermussmann/Documents/GitHub/02427-AdvancedTimeSeries-CEX1/code_matlab';
copyfile(current_file, new_directory);