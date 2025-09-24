%% Value Function Iteration (Stochastic) %%
clear all; close all; clc
% parameters
beta = 0.95;    % discount factor
alpha = 0.33;   % capital share
delta = 0.1;    % depreciation rate
gamma = 2;      % Inverse of EIS
kappa = 0.8;

rho = 0.95;     % persistence of shock
nZ = 7;         % # of shock nodes
mu = 0;         % unconditional mean of process
m = 3;          % cover
sigma = 0.01;   % std of shock

nK = 200;       % # of capital grid
tol = 1e-7;     % tolerance
iter_max = 1000;      
err_V = 999;

% steady-state
Kss = ((1 - beta*(1 - delta))/(alpha*beta))^(1/(alpha - 1)) ;
Css = Kss^alpha - delta*Kss ;

% Discretize Shocks
[Z, Zprob] = tauchen(nZ, mu, rho, sigma, m);
Zgrid = exp(Z);

% Step 1: Set the grid space for K
Kmin = (1-kappa)*Kss;
Kmax = (1+kappa)*Kss;
Kgrid = linspace(Kmin, Kmax, nK)';

% Step 2: Initialization
V = zeros(nK,nZ);     % Value function
V0 = zeros(nK,nZ);    % Old value function
Kp = zeros(nK,nZ);    % policy function (k_prime)
Kpi = zeros(nK,nZ);   % index of k_prime
Cons = zeros(nK,nZ);  % policy function (consumption)
W = zeros(nK,nZ);   

tic
%% Step 3: Iterate Value Function until convergence
for iter = 1:iter_max
    W = beta*V0*Zprob';
    for iZ = 1:nZ
        % Monotonicity
        iKp0 = 1;
        for iK = 1:nK
            V(iK,iZ) = -1e7;
            for iKp = iKp0:nK
                c_val = Zgrid(iZ)*Kgrid(iK)^alpha + (1 - delta)*Kgrid(iK) - Kgrid(iKp);
                if c_val < 0
                    u_val = -10000;
                else
                    u_val = c_val^(1 - gamma)/(1 - gamma);
                end
                
                v_val = u_val + W(iKp,iZ);
                
                if v_val > V(iK,iZ)
                    V(iK,iZ) = v_val;
                    Cons(iK,iZ) = c_val;
                    Kp(iK,iZ) = Kgrid(iKp);
                    iKp0 = iKp;
                    Kpi(iK,iZ) = iKp;
                % Concavity
                elseif v_val < V(iK,iZ)
                    break
                end               
            end           
        end       
    end
    
    err_V = max(max(abs(V - V0)))/max(max(abs(V0)));
    fprintf('iter: %d     error:%1.8f\n',iter,err_V)
    if err_V<tol
        break
    end
    
    % update
    V0 = V;    
    
end

toc



