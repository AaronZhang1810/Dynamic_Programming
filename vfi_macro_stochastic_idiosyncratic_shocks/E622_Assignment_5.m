%% Author: Jihuan Zhang 'Aaron'

%% Parameters Setup
% Housekeeping
clear all; close all; clc

alpha = 1/3; % capital share
beta = 0.96; % discount factor
delta = 0.065; % depreciation
sigma = 2; % parameter for the utility function

% Capital grid setup
na = 200; % number of capital nodes
kss = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1)); % steady-state capital
a_max = 5*kss; % upper bound for capital grid
b = 0; % customized lower bound by model specification
agrid = linspace(b,a_max,na);

% Idiosyncratic shock grid setup
mu = 0; % parameter for the distribution of aggregate shocks
sigma_v = 0.2; % parameter for the distribution of aggregate shocks
rho = 0.9; % persistency of aggregate shocks
ne = 5; % number of shock nodes
m = 3; % max +- std. devs.
% simulating a Markov chain that approximates the AR(1) process
[egrid, eprob] = tauchen(ne, mu, rho, sigma_v, m);
% the transition matrix specifies transition probability from row to column
egrid = exp(egrid); % since we assume log(e) ~ AR(1)

I_max = 1000; % number of households to be simulated
T_max = 2000; % number of periods to be simulated
e_index = zeros(I_max,T_max); % initialize state indices
e = zeros(I_max,T_max); % initialize idiosyncratic shocks
mc = dtmc(eprob); % create a discrete-time Markov chain
rng('default'); % reset random seed
for i = 1:I_max % simulation
    rng(i); % fix seed for each iteration for reproducing results
    e_val = simulate(mc,T_max-1); % simulate an individual for T periods
    e_index(i,:) = e_val; % store the simulated state indices
    e(i,:) = egrid(e_val); % store the simulated idiosyncratic shocks
end

% Iteration setup - outer loop
tol_r = 1e-3;
distance_r = 1e3;
r0 = 1/beta-1-1e-5; % initial guess of r around the deterministic steady state
% according to the theory, optimal r is in (-delta, 1/beta-1)
r1 = 1/beta-1-1e-5;

% Iteration setup - inner loop
V0 = zeros(na,ne); % old value function
V1 = zeros(na,ne); % new value function
a_policy = zeros(na,ne); % policy function (savings)
a_index = zeros(na,ne); % index of policy function
c_policy = zeros(na,ne); % policy function (consumption)
v = zeros(na,1);
itermax = 1000;
tol_v = 1e-3;

%% Solve the model - standard algorithm
% this algorithm is standard and easily understandable but takes time
% we can also make use of the monotonicity of policy functions and
% concavity of the value function to improve the algorithm, as I did for
% a similar problem with aggregate shocks
while distance_r > tol_r
    r0 = 0.1*r1 + 0.9*r0; % slow updating rule
    w0 = (1-alpha)*(alpha/(r0+delta))^(alpha/(1-alpha)); % solve for wage
    distance_v = 1e3; % the distance should be reset after the completion of each vfi
    % vfi for each updating of r
    for iter = 1:itermax
        for ie = 1:ne
            for ia = 1:na
                for iap = 1:na
                    c = w0*egrid(ie)+(1+r0)*agrid(ia)-agrid(iap);
                    ev = beta*sum(eprob(ie,:).*V0(iap,:));
                    if c <= 0 % penalize unrealistic consumption
                        v(iap) = -1e5;
                    else
                        if sigma == 1 % when sigma=1, utility function is logarithm
                            v(iap) = log(c)+ev;
                        else
                            v(iap) = (c^(1-sigma)-1)/(1-sigma)+ev;
                        end
                    end
                end
                [V1(ia,ie), a_index(ia,ie)] = max(v);
            end
        end
        % we compute the vfi error by the relevant maximum distance of all
        % entries. Note that "max(abs(V1-V0))" gives a row vector, we need
        % double max
        distance_v = max(max(abs(V1-V0)))/max(max(abs(V0)));
        fprintf('iter: %d error:%1.8f\n current r: %1.5f\n',iter,distance_v,r0);
        disp("distance:" + distance_r);
        if distance_v < tol_v
            break
        end
        V0 = V1; % updating
    end

    % Recover the policy functions
    for ia = 1:na
        for ie = 1:ne
            a_policy(ia,ie) = agrid(a_index(ia,ie));
            c_policy(ia,ie) = w0*egrid(ie)+(1+r0)*agrid(ia)-a_policy(ia,ie);
        end
    end

    % Individual Simulation
    % we just solved for the optimal functions
    % note that the functional forms are the same for every HH at every period
    % now we compute the realized decision and utility paths for each HH
    % we accomplish this using the simulated idiosyncratic shocks
    % at = zeros(I,T);
    I = 1000;
    T = 500; % we use the first 500 periods
    et = e(1:I,1:T); eit = e_index(1:I,1:T); % we use a subsample of the Markov chain
    at = ones(I,T+1)*(-b); % T+1 includes zero as the starting time
    ait = ones(I,T+1); % assume every HH starts with the lowest possible endowment
    values = zeros(I,T);
    for i = 1:I
        for t = 1:T
            ei = eit(i,t); % the index of shock state for household i at time t
            ai = ait(i,t); % the index of capital state for household i at time t
            % a_index(i,j) means that if a=a(grid=i) and e=e(grid=j), the 
            % best response of household i at time t is to choose a'=a(grid=a_index(i,j))
            api = a_index(ai,ei); % asset decision at t according to the policy function
            aprime = a_policy(ai,ei); % the desicion of HH i at time t
            values(i,t) = V1(ai,ei); % the optimal expected utility of HH i at time t
            % updating
        	ait(i,t+1) = api; % next period capital state index is the optimal
            % capital index of current period
            at(i,t+1) = aprime; % next period capital state is the optimal
            % capital decision of current period
        end
    end
    % calculate moments, i.e., average capital and labor at the steady state
    K = mean(mean(at(:,T-10:T))); % use the last 10 periods because they are 
    % closer to the steady-state values
    N = mean(mean(e(:,T-10:T))); % steady-state average labor

    % then use the moments to update r1, using firms' FOC
    r1 = alpha*(K/N)^(alpha-1)-delta; % calculating mean instead of aggregates
    % won't hurt because the population will always be cancelled
    distance_r = abs(r0-r1); % one may also consider relevant distance
end

%% Ploting
legend_str = ["lowest","median","highest"];
figure(1)
hold on
plot(agrid,a_policy(:,1));
plot(agrid,a_policy(:,3));
plot(agrid,a_policy(:,5));
ylabel 'Policy Function for Savings';
legend(legend_str,'Location','southeast');
hold off
print(gcf,'E622_Assignment_5_Policy_Function_Savings.png','-dpng','-r600');

figure(2)
hold on
plot(agrid,c_policy(:,1));
plot(agrid,c_policy(:,3));
plot(agrid,c_policy(:,5));
ylabel 'Policy Function for Consumptions';
legend(legend_str,'Location','southeast');
hold off
print(gcf,'E622_Assignment_5_Policy_Function_Consumptions.png','-dpng','-r600');

figure(3)
hold on
plot(agrid,V1(:,1));
plot(agrid,V1(:,3));
plot(agrid,V1(:,5));
ylabel 'Value Function';
legend(legend_str,'Location','southeast');
hold off
print(gcf,'E622_Assignment_5_Value_Function.png','-dpng','-r600');