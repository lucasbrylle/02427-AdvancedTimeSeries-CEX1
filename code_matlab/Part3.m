clc; clear; close all;

% OBS!!!!
% Integral of the theoretical conditional means (Eq. 3-4 in book)
% Theoretical conditional mean of the SETAR(2,1,1) is calculated in Part 2


% Script for calculation of a confidence bands using the cumulative means technique.
% First generate x: a time series with a realization of a Markov process, i.e. some
% model where x_t is dependent on x_{t-1}


% ---------------- Change in threshold ------------------------------------

% Parameters
T = 3000; % Number of observations
threshold = [-1, 1]; % Threshold for the SETAR model
epsilon = randn(T, 1); % Noise term
seed = 5;
phi = 0.7; % AR coefficient for both regimes
c1 = 4; % Constant for regime 1
c2 = -4; % Constant for regime 2
sigma2 = 1; % Variance of white noise


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
    n_bin = 50;
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


    % Calculate Theoretical Cumulative Mean using bin points
    theoretical_cum_mean = zeros(n_bin,1);
    
    for i = 1:n_bin
        if bin_points(i) <= threshold(j)
            % Regime 1 cumulative mean
            theoretical_cum_mean(i) = c1 * (bin_points(i) - bin(1)) + (phi / 2) * (bin_points(i)^2 - bin(1)^2);
        else
            % Regime 2 cumulative mean
            theoretical_cum_mean(i) = c1 * (threshold(j) - bin(1)) + ...
                                      (phi / 2) * (threshold(j)^2 - bin(1)^2) + ...
                                      c2 * (bin_points(i) - threshold(j)) + ...
                                      (phi / 2) * (bin_points(i)^2 - threshold(j)^2);
        end
    end
    
    % Save theoretical cumulative mean for plotting
    theoretical_cum_mean_struc{j} = theoretical_cum_mean;


end

% Plotting

% Define the custom color
customRed = [195/255, 31/255, 51/255];
customBlue = [24/255, 54/255, 104/255];
customGreen = [60/255, 179/255, 113/255];

% Font size and weight
fontSize = 13;
fontWeight = 'bold';

figure('Position',[100, 100, 800, 800])
subplot(2,1,1)
plot(bin_points, Lambda_struc{1}, "color", customBlue, "linewidth", 2);
hold on
plot(bin_points, Lambda_low_struc{1}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points, Lambda_up_struc{1}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points, theoretical_cum_mean_struc{1}, 'color', customGreen, 'linewidth', 2);
hold off
title('Threshold = -1', 'FontWeight', fontWeight, 'FontSize',fontSize)
ylabel('$\mathbf{\Lambda}$', 'Interpreter', 'latex');
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;


subplot(2,1,2)
plot(bin_points, Lambda_struc{2}, "color", customBlue, "linewidth", 2);
hold on
plot(bin_points, Lambda_low_struc{2}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points, Lambda_up_struc{2}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points, theoretical_cum_mean_struc{2}, 'color', customGreen, 'linewidth', 2);
hold off
ylabel('$\mathbf{\Lambda}$', 'Interpreter', 'latex');
xlabel('Bin', 'FontWeight', fontWeight, 'FontSize',fontSize);
title('Threshold = 1', 'FontWeight', fontWeight, 'FontSize',fontSize)
legend('Estimated Mean', 'Lower Bound', 'Upper Bound', 'Theoretical Mean', 'Location','northwest')
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;


% ------------------ Change in number of bins -----------------------------

% Parameters
T = 3000; % Number of observations
threshold = 0; % Threshold for the SETAR model
epsilon = randn(T, 1); % Noise term
seed = 5;
phi = 0.7; % AR coefficient for both regimes
c1 = 4; % Constant for regime 1
c2 = -4; % Constant for regime 2
sigma2 = 1; % Variance of white noise

%%--------
% Hist regression prm.
bin = [-2 2];
n_bin = [10 50];
%%--------


for j = 1:length(n_bin)
    % Initialization
    x = zeros(T, 1); % Initialize the time series
    x(1) = epsilon(1);
    % SETAR(2,1,1) Model Simulation
    for t = 2:T
        if x(t-1) <= threshold
            x(t) = 4 + 0.7 * x(t-1) + epsilon(t);       % Regime 1
        else
            x(t) = -4 + 0.7 * x(t-1) + epsilon(t);      % Regime 2
        end
    end
        
    %%--------
    % Hist regression prm.
    h = (bin(2) - bin(1))/n_bin(j);
    bin_points = (bin(1)+h/2):h:(bin(2)-h/2);
    %%--------


    
    %%--------
    % Hist regression
    % the starting bin
    cur_bin=[bin_points(1)-0.5*h bin_points(1)+0.5*h];
    % init
    lambda=zeros(n_bin(j),1);
    gamma=zeros(n_bin(j),1);
    f_hat=zeros(n_bin(j),1);
    
    for i=1:n_bin(j)
        index=(x(1:end-1)>cur_bin(1) & x(1:end-1)<=cur_bin(2));
        if (sum(index)>5)
            lambda(i) = sum( x(2:end).*index ) / sum(index);
            % f_hat is needed for confidence band calculation
            f_hat(i) = (n_bin(j)*h)^(-1) * sum(index);
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
    h_hat = zeros(n_bin(j),1);
    for i=1:n_bin(j)
        h_hat(i) = gamma(i)/f_hat(i);
    end
    H_hat = cumsum(h_hat*h);
    
    H_hat_b = H_hat(n_bin(j));
    Lambda_lower = Lambda - c .* n_bin(j).^(-0.5) .* H_hat_b.^(0.5) .* (1 + H_hat/H_hat_b);
    Lambda_upper = Lambda + c .* n_bin(j).^(-0.5) .* H_hat_b.^(0.5) .* (1 + H_hat/H_hat_b);
    %%--------
    % Save data
    Lambda_struc{j} = Lambda;
    Lambda_low_struc{j} = Lambda_lower;
    Lambda_up_struc{j} = Lambda_upper;


    % Calculate Theoretical Cumulative Mean using bin points
    theoretical_cum_mean = zeros(n_bin(j),1);
    
    for i = 1:n_bin(j)
        if bin_points(i) <= threshold
            % Regime 1 cumulative mean
            theoretical_cum_mean(i) = c1 * (bin_points(i) - bin(1)) + (phi / 2) * (bin_points(i)^2 - bin(1)^2);
        else
            % Regime 2 cumulative mean
            theoretical_cum_mean(i) = c1 * (threshold - bin(1)) + ...
                                      (phi / 2) * (threshold^2 - bin(1)^2) + ...
                                      c2 * (bin_points(i) - threshold) + ...
                                      (phi / 2) * (bin_points(i)^2 - threshold^2);
        end
    end
    
    % Save theoretical cumulative mean for plotting
    theoretical_cum_mean_struc{j} = theoretical_cum_mean;

    % Save binpoints
    bin_points_struc{j} = bin_points;

end

% Plotting

% Define the custom color
customRed = [195/255, 31/255, 51/255];
customBlue = [24/255, 54/255, 104/255];
customGreen = [60/255, 179/255, 113/255];

% Font size and weight
fontSize = 13;
fontWeight = 'bold';

figure('Position',[100, 100, 800, 800])
subplot(2,1,1)
plot(bin_points_struc{1}, Lambda_struc{1}, "color", customBlue, "linewidth", 2);
hold on
plot(bin_points_struc{1}, Lambda_low_struc{1}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points_struc{1}, Lambda_up_struc{1}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points_struc{1}, theoretical_cum_mean_struc{1}, 'color', customGreen, 'linewidth', 2);
hold off
title('Number of bins = 10', 'FontWeight', fontWeight, 'FontSize',fontSize)
ylabel('$\mathbf{\Lambda}$', 'Interpreter', 'latex');
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;


subplot(2,1,2)
plot(bin_points_struc{2}, Lambda_struc{2}, "color", customBlue, "linewidth", 2);
hold on
plot(bin_points_struc{2}, Lambda_low_struc{2}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points_struc{2}, Lambda_up_struc{2}, "color", customRed, "linewidth", 1, "LineStyle","--");
hold on
plot(bin_points_struc{2}, theoretical_cum_mean_struc{2}, 'color', customGreen, 'linewidth', 2);
hold off
ylabel('$\mathbf{\Lambda}$', 'Interpreter', 'latex');
xlabel('Bin', 'FontWeight', fontWeight, 'FontSize',fontSize);
title('Number of bins = 50', 'FontWeight', fontWeight, 'FontSize',fontSize)
legend('Estimated Mean', 'Lower Bound', 'Upper Bound', 'Theoretical Mean', 'Location','northwest')
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

current_file = 'Part3.m';
new_directory = '/home/olivermussmann/Documents/GitHub/02427-AdvancedTimeSeries-CEX1/code_matlab';
copyfile(current_file, new_directory);