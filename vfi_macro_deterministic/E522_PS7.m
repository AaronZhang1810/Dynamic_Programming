%% Author: Jihuan Zhang 'Aaron'

%% Settings
close all; clear all; clc;

% The model (representative agent, no leisure and nonstochastic)
% utility function: u = log(c)
% production function: F(K,1) = K^alpha * 1^(1-alpha)
% aggregate resource constraint: c + kprime = k^alpha + (1-delta)*k
% steady-state capital: kss = ((1/beta -1+delta)/alpha)^(1/(alpha-1))

% Parameters Setup
alpha = 0.33; % share of capital
beta = 0.96; % discount factor
delta = 1; % depreciation rate of capital

% Capital grid setup
kappa = 0.9; % the bound for capital grid
kss = ((1/beta -1+delta)/alpha)^(1/(alpha-1)); % steady-state capital
na = 200; % the number of grid points
k_grid = linspace((1-kappa)*kss,(1+kappa)*kss,na); % capital grid (row vector)

% Value function and policy function setup
v_old = zeros(na,1); % create a column vector for V0. Will be updated later.
v_new = zeros(na,1); % create a column vector for V1. Will be updated later.
v_val = zeros(na,1); % create a column vector to store value to find maximizers.
k_policy = zeros(na,1); % create a column vector for policy function for capital.
c_policy = zeros(na,1); % create a column vector for policy function for consumption.
optimal_position = zeros(na,1); % create a column vector to store index of maximizers.

% Iteration settings
tolerance = 1e-5; % tolerance level
iter_max = 2000; % maximum number of iteration
distance = tolerance+1;

%% Value function iteration
for iter = 1:iter_max
    for ik = 1:na
    	for ikprime = 1:na
        	consumption = k_grid(ik)^alpha + (1-delta)*k_grid(ik) - k_grid(ikprime);
        	if consumption > 0
        	    v_val(ikprime) = log(consumption) + beta*v_old(ikprime,1);
            else
            	v_val(ikprime) = -100000;
            end
        end
        [v_new(ik,1), optimal_position(ik,1)] = max(v_val);
    end
    distance = max(abs((v_new-v_old)./v_old)); % this is the relative distance which is more precise
    v_old = v_new;
    disp("Current distance:" + distance)
    if distance <= tolerance
        break
    end
end

% Policy function
for i = 1:na
    k_policy(i,1) = k_grid(optimal_position(i,1)); % k prime is a function of k
    c_policy(i,1) = k_grid(i)^alpha + (1-delta)*k_grid(i) - k_policy(i,1); % c is a function of k
end

%% Plot
figure(1)
plot(k_grid,v_new);
ylabel 'Value Function';

figure(2)
plot(k_grid,k_policy);
ylabel('Policy function for capital');

figure(3)
plot(k_grid,c_policy);
ylabel('Policy function for consumption');