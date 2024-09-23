clc; clear; close all;

% OBS!!!!
% Integral of the theoretical conditional means (Eq. 3-4 in book)
% Theoretical conditional mean of the SETAR(2,1,1) is calculated in Part 2


% Script for calculation of a confidence bands using the cumulative means technique.
% First generate x: a time series with a realization of a Markov process, i.e. some
% model where x_t is dependent on x_{t-1}

% Parameters
T = 1000; % Number of observations
threshold = [0, 1, 2]; % Threshold for the SETAR model
epsilon = randn(T, 1); % Noise term
seed = 5;


%%

for j = 1:length(threshold)
    % Initialization
    x = zeros(T, 1); % Initialize the time series
    x(1) = epsilon(1);
    % SETAR(2,1,1) Model Simulation
    for t = 2:T
        if x(t-1) <= threshold(j)
            x(t) = 4 + 0.7 * x(t-1) + epsilon(t);       % Regime 1
        else
            x(t) = -4 + 0.7 * x(t-1) + epsilon(t);      % Regime 2
        end
    end
    
    %%--------
    % Hist regression prm.
    bin = [-2 2];
    n_bin = 20;
    h = (bin(2) - bin(1))/n_bin;
    bin_points = (bin(1)+h/2):h:(bin(2)-h/2);
    %%--------
    
    
    %%--------
    % Hist regression
    % the starting bin
    cur_bin=[bin_points(1)-0.5*h bin_points(1)+0.5*h];
    % init
    lambda=zeros(n_bin,1);
    gamma=zeros(n_bin,1);
    f_hat=zeros(n_bin,1);
    
    for i=1:n_bin
        index=(x(1:end-1)>cur_bin(1) & x(1:end-1)<=cur_bin(2));
        if (sum(index)>5)
            lambda(i) = sum( x(2:end).*index ) / sum(index);
            % f_hat is needed for confidence band calculation
            f_hat(i) = (n_bin*h)^(-1) * sum(index);
            gamma(i) = sum( (x(2:end) - lambda(i) ).^2 .* index ) / sum(index);
        else
            fprintf('Need more points')
            break
        end
    
        % move to next bin
        cur_bin=cur_bin+h;
    end
    
    
    % Make confidence bands for the cumulated function. Def. (3.10).
    % 95% confidence band, c is found in table 3.1
    c = 1.273;
    
    Lambda = cumsum(lambda*h);
    h_hat = zeros(n_bin,1);
    for i=1:n_bin
        h_hat(i) = gamma(i)/f_hat(i);
    end
    H_hat = cumsum(h_hat*h);
    
    H_hat_b = H_hat(n_bin);
    Lambda_lower = Lambda - c .* n_bin.^(-0.5) .* H_hat_b.^(0.5) .* (1 + H_hat/H_hat_b);
    Lambda_upper = Lambda + c .* n_bin.^(-0.5) .* H_hat_b.^(0.5) .* (1 + H_hat/H_hat_b);
    %%--------
    % Save data
    Lambda_struc{j} = Lambda;
    Lambda_low_struc{j} = Lambda_lower;
    Lambda_up_struc{j} = Lambda_upper;


end

% Plotting

% Define the custom color
customRed = [195/255, 31/255, 51/255];
customBlue = [24/255, 54/255, 104/255];

% Font size and weight
fontSize = 13;
fontWeight = 'bold';

figure
subplot(3,1,1)
plot(bin_points, Lambda_struc{1}, "color", customBlue, "linewidth", 2);
hold on
plot(bin_points, Lambda_low_struc{1}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points, Lambda_up_struc{1}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold off
xlabel('Bin', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('CMean', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;


subplot(3,1,2)
plot(bin_points, Lambda_struc{2}, "color", customBlue, "linewidth", 2);
hold on
plot(bin_points, Lambda_low_struc{2}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points, Lambda_up_struc{2}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold off
xlabel('Bin', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('CMean', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

subplot(3,1,3)
plot(bin_points, Lambda_struc{3}, "color", customBlue, "linewidth", 2);
hold on
plot(bin_points, Lambda_low_struc{3}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points, Lambda_up_struc{3}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold off
ylabel('$\mathbf{\hat{\Lambda}}$, $\mathbf{\Lambda}$ ', 'Interpreter', 'latex');
xlabel('Bin', 'FontWeight', fontWeight, 'FontSize', fontSize);
%ylabel('CMean', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

current_file = 'Part3.m';
new_directory = '/home/olivermussmann/Documents/GitHub/02427-AdvancedTimeSeries-CEX1/code_matlab';
copyfile(current_file, new_directory);