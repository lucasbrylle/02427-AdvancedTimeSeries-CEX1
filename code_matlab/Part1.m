clc; clear; close all;

% Simulating SETAR(2,1;1,2), IGAR(2;1,2), and MMAR(2;1) Models
T = 2000; % Number of observations

% Initialize time series
x_setar = zeros(T, 1);
x_igar = zeros(T, 1);
x_mmar = zeros(T, 1);

% Noise term
epsilon = randn(T, 1);

% ------------------- SETAR(2,1;1,2) model --------------------------------

% Threshold
threshold = 0.5;

% Counter
Reg_SETAR = zeros(T, 1);

% Simulation
for t = 2:T
    if x_setar(t-1) <= threshold
        x_setar(t) = 4 + 0.7 * x_setar(t-1) + epsilon(t);
        Reg_SETAR(t, 1) = 1;
    else
        x_setar(t) = -4 + 0.7 * x_setar(t-1) + epsilon(t);
        Reg_SETAR(t, 1) = 2;
    end
end

% Display Regime 1 percentage
disp(['Percentage of points in Regime 1: ', num2str(sum(Reg_SETAR==1)/length(Reg_SETAR) * 100), '%']);


% --------------------- IGAR(2;1,2) model ---------------------------------

% Counter
Reg_IGAR = zeros(T, 1);

% Probability for Regime 1
prob = 0.75;

% Simulation
for t = 2:T
    if rand < prob
        x_igar(t) = 4 + 0.7 * x_igar(t-1) + epsilon(t);
        Reg_IGAR(t, 1) = 1;
    else
        x_igar(t) = -4 + 0.7 * x_igar(t-1) + epsilon(t);
        Reg_IGAR(t, 1) = 2;
    end
end

% Display Regime 1 percentage
disp(['Percentage of points in Regime 1: ', num2str(sum(Reg_IGAR==1)/length(Reg_IGAR) * 100), '%']);


% ------- MMAR(2;1) model with Markov chain regime switching --------------

% Counter
state_count = zeros(T, 1);


% Transition matrix
P = [0.6, 0.4; 0.1, 0.9];

% P(1, 1) = 0.6: The probability of staying in state 1 if you are currently in state 1.
% P(1, 2) = 0.4: The probability of transitioning to state 2 if you are currently in state 1.
% P(2, 1) = 0.1: The probability of transitioning to state 1 if you are currently in state 2.
% P(2, 2) = 0.9: The probability of staying in state 2 if you are currently in state 2.

% Solve for the stationary distribution
pi = null(P' - eye(size(P)),'r'); % 'r' ensures the solution is real
pi = pi / sum(pi); % Normalize to ensure it sums to 1

% Display the stationary distribution
disp('Stationary distribution:');
disp(pi');

% Initialization
state = 1;

% Simulation
for t = 2:T
    if state == 1
        x_mmar(t) = 4 + 0.7 * x_mmar(t-1) + epsilon(t);
    else
        x_mmar(t) = -4 + 0.7 * x_mmar(t-1) + epsilon(t);
    end
    % Transition to next state
    state = find(rand < cumsum(P(state, :)), 1);
    state_count(t, 1) = state;
end

% Display Regime 1 percentage
disp(['Percentage of points in Regime 1: ', num2str(sum(state_count==1)/length(state_count) * 100), '%']);


%% ------------------------- Plotting -------------------------------------

% ------------------------ Initilization ----------------------------------

% Time vector
time = (1:T)';

% Define the custom color
customRed = [195/255, 31/255, 51/255];
customBlue = [24/255, 54/255, 104/255];

% Font size and weight
fontSize = 13;
fontWeight = 'bold';

% -------------------- Plotting the simulations ---------------------------

% Plot SETAR model with small custom-colored scatters
figure;
subplot(2, 1, 1);
scatter(time(1:1000), x_setar(1:1000), 10, customBlue); % Small custom-colored scatter points
xlabel('t', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('X_t', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

subplot(2, 1, 2);
scatter(x_setar(1:end-1), x_setar(2:end), 10, customRed); % Small custom-colored scatter points
xlabel('X_{t-1}', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('X_t', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

% Save the plot as an SVG file
saveas(gcf, 'Task1_SETAR.svg');

% Counter plot
figure;
scatter(time(1:50), Reg_SETAR(1:50), 10, customBlue, "filled"); % Small custom-colored scatter points
xlabel('Index', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('Regime', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
ylim([0.95 2.05])
box on;

% Save the plot as an SVG file
saveas(gcf, 'Task1_SETAR_count.svg');


% % Plot IGAR model
figure;
subplot(2, 1, 1);
scatter(time(1:1000), x_igar(1:1000), 10, customBlue); % Small custom-colored scatter points
xlabel('t', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('X_t', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

subplot(2, 1, 2);
scatter(x_igar(1:end-1), x_igar(2:end), 10, customRed); % Small custom-colored scatter points
xlabel('X_{t-1}', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('X_t', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

% Save the plot as an SVG file
saveas(gcf, 'Task1_IGAR.svg');

% % Counter plot
figure;
scatter(time(1:50), Reg_IGAR(1:50), 10, customBlue, "filled"); % Small custom-colored scatter points
xlabel('Index', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('Regime', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
ylim([0.95 2.05])
box on;

% Save the plot as an SVG file
saveas(gcf, 'Task1_IGAR_count.svg');
% 
% % 
% % % Plot MMAR model
figure;
subplot(2, 1, 1);
scatter(time(1:1000), x_mmar(1:1000), 10, customBlue); % Small custom-colored scatter points
xlabel('t', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('X_t', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

subplot(2, 1, 2);
scatter(x_mmar(1:end-1), x_mmar(2:end), 10, customRed); % Small custom-colored scatter points
xlabel('X_{t-1}', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('X_t', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
box on;

% Save the plot as an SVG file
saveas(gcf, 'Task1_MMAR.svg');

% % Counter plot
figure;
scatter(time(1:50), state_count(1:50), 10, customBlue, "filled"); % Small custom-colored scatter points
xlabel('Index', 'FontWeight', fontWeight, 'FontSize', fontSize);
ylabel('Regime', 'FontWeight', fontWeight, 'FontSize', fontSize);
set(gca, 'FontWeight', fontWeight, 'FontSize', fontSize); % Set font size and weight for ticks
ylim([0.95 2.05])
box on;

% Save the plot as an SVG file
saveas(gcf, 'Task1_MMAR_count.svg');

current_file = 'Part1.m';
new_directory = '/home/olivermussmann/Documents/GitHub/02427-AdvancedTimeSeries-CEX1/code_matlab';
copyfile(current_file, new_directory);