%% Author: Jihuan Zhang 'Aaron'

%% Parameters Setup
% Housekeeping
clear all; close all; clc

alpha = 1/3; % capital share
beta = 143/144; % discount factor
delta = 1/48; % depreciation
sigma = 2; % parameter for the utility function

% Capital grid setup
kappa = 0.4;
nk = 200; % number of capital nodes
kss = (alpha/((1/beta)+delta-1))^(1/(1-alpha)); % steady state capital
kgrid = linspace((1-kappa)*kss, (1+kappa)*kss, nk);

% Aggregate shock grid setup
mu = 0; % parameter for the distribution of aggregate shocks
sigma_e = 0.007; % parameter for the distribution of aggregate shocks
rho = 0.975; % persistency of aggregate shocks
nz = 7; % number of shock nodes
m = 3; % max +- std. devs.
[zgrid, zprob] = tauchen(nz, mu, rho, sigma_e, m); % simulating Markov chain
% that approximates the AR(1) process
% the transition matrix specifies transition probability from row to column
zgrid = exp(zgrid); % since we assume log(z) ~ AR(1)

% Iteration setup
tol = 1e-7; % tolerance
iter_max = 1000;
V0 = zeros(nk,nz);
V1 = zeros(nk,nz);
Vnew = zeros(nk,1);
C_policy = zeros(nk,nz); % consumption policy function
K_policy = zeros(nk,nz); % capital policy function
optimal_index = zeros(nk,nz); % optimal grid position

%% Value function iteration - Standard Algorithm

% this algorithm is standard and easily understandable but takes time
% the next block is a faster algorithm but restricted to a smaller subset
% of dynamic programming problems.
for iter = 1:iter_max
    for iz = 1:nz
        for ik = 1:nk
            for ikp = 1:nk
                c  = zgrid(iz)*kgrid(ik)^alpha+(1-delta)*kgrid(ik)-kgrid(ikp);
                if c > 0
                    Vnew(ikp) = (c^(1-sigma)-1)/(1-sigma) + beta*sum(zprob(iz,:).*V0(ikp,:));
                else
                    Vnew(ikp) = -100000;
                end
            end
            [V1(ik,iz), optimal_index(ik,iz)] = max(Vnew);
        end
    end
    % we compute the vfi error by the relevant maximum distance of all
    % entries. Note that "max(abs(V1-V0))" gives a row vector, we need
    % double max
    distance = max(max(abs((V1-V0)))/max(max(abs(V0))));
    fprintf('iter: %d error:%1.8f\n',iter,distance);
    if distance < tol
        break
    end
    % updating
    V0 = V1;
end

% Recover the policy functions
for ik = 1:nk
    for iz = 1:nz
        indx = optimal_index(ik,iz);
        K_policy(ik,iz) = kgrid(indx);
        C_policy(ik,iz) = zgrid(iz)*kgrid(ik)^alpha+(1-delta)*kgrid(ik)-kgrid(indx);
    end
end

%% Value function iteration - Faster Algorithm

% this algorithm works only for cases where policy functions are nondecreasing
% and value function is concave.
tic
for iter = 1:iter_max
    W = beta*V0*zprob'; % W(i,j) is the expectation of beta*V(z',k') at
    % k' = i and z = j. We will later find maximum for each column, i.e., z
    for iz = 1:nz
        % Monotonicity: policy functions are increasing in k
        % e.g., if k'(grid=10) is the maximizer for k(grid=1),
        % then k'(grid=9) cannot be the maximizer for k(grid=2)
        % and this is true for any fixed state z
        % using this observation, we can significantly imporve the speed
        ikp0 = 1;
        for ik = 1:nk
            V1(ik,iz) = -1e7;
            for ikp = ikp0:nk
                c = zgrid(iz)*kgrid(ik)^alpha+(1-delta)*kgrid(ik)-kgrid(ikp);
                if c < 0
                    v = -1e-5;
                else
                    % modification is needed if sigma is set to be one
                    v = (c^(1-sigma)-1)/(1-sigma) + W(ikp,iz);
                end
                if v > V1(ik,iz)
                    V1(ik,iz) = v;
                    C_policy(ik,iz) = c;
                    K_policy(ik,iz) = kgrid(ikp);
                    ikp0 = ikp; % so that for next k, we start from iKp
                    optimal_index(ik,iz) = ikp;
                % Concavity: value function increases with k' first but
                % then declines with k', otherwise we would always pick k'
                % as large as possible.
                elseif v < V1(ik,iz)
                    % if v starts decreasing, then the previous v is the maximum
                    break % we return immediately back to the next state k
                    % for the next k, we don't start at grid=1, instead, we
                    % start at the last optimal capital node, by
                    % monotonicity
                end
            end
        end
    end
    % we compute the vfi error by the relevant maximum distance of all
    % entries. Note that "max(abs(V1-V0))" gives a row vector, we need
    % double max
    distance = max(max(abs(V1-V0)))/max(max(abs(V0)));
    fprintf('iter: %d error:%1.8f\n',iter,distance);
    if distance < tol
        break
    end
    % updating
    V0 = V1;
end
toc

%% Plot Value Function
legend_str = string(zeros(nz,1));
figure(1)
hold on
for i = 1:nz
    plot(kgrid,V0(:,i));
    legend_str{i} = ['z = ' num2str(zgrid(i),4)];
end
xlabel 'Nk = 200, kappa = 0.4';
ylabel 'Value Function';
legend(legend_str,'Location','southeast');
hold off
print(gcf,'E622_Assignment_3_Value_Function_nk=200_kappa=0.4.png','-dpng','-r600');
%% Plot Policy Functions
figure(2)
plot(kgrid,C_policy);
xlabel 'Nk = 200, kappa = 0.4';
ylabel 'Next Period Consumption';
legend(legend_str,'Location','southeast');
print(gcf,'E622_Assignment_3_Policy_Functions_Consumption_nk=200_kappa=0.4.png','-dpng','-r600');

figure(3)
plot(kgrid,K_policy);
xlabel 'Nk = 200, kappa = 0.4';
ylabel 'Next Period Capital';
legend(legend_str,'Location','southeast');
print(gcf,'E622_Assignment_3_Policy_Functions_Capital_nk=200_kappa=0.4.png','-dpng','-r600');
hold off