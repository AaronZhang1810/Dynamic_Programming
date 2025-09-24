%% State space representation
% S(t) = A*S(t-1) + B*e(t)
% X(t) = C*S(t-1) + D*e(t)

clc
clear 
close all

dynare model1

p_c = 1;
p_k = 2;
p_a = 3;

A =[oo_.dr.ghx(oo_.dr.inv_order_var(p_k),:)
    oo_.dr.ghx(oo_.dr.inv_order_var(p_a),:)];

B =[oo_.dr.ghu(oo_.dr.inv_order_var(p_k),:)
    oo_.dr.ghu(oo_.dr.inv_order_var(p_a),:)];

C = oo_.dr.ghx(oo_.dr.inv_order_var(p_c),:);

D = oo_.dr.ghu(oo_.dr.inv_order_var(p_c),:);






%% Value Function Iteration (Stochastic) %%

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

% Discreteize Shocks
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
% Step 3: Iterate Value Function until convergence
for iter = 1:iter_max
    W = beta*V0*Zprob';
    for iK = 1:nK
        % Monotonicity
        iKp0 = 1;
        if iK>1 
            iKp0 = max(Kpi(iK-1,iZ),1);
        end
        for iZ = 1:nZ
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

cons = Css*(C(1,1)*(Kgrid-Kss)/Kss + D(1,1)*(Zgrid(3)-1)) + Css;
kprime = Kss*(A(1,1)*(Kgrid-Kss)/Kss + B(1,1)*(Zgrid(3)-1)) + Kss;

kprime_VFI = Kp(:,3);
cons_VFI = Cons(:,3);


figure()
plot(Kgrid,cons,'LineWidth',2)
hold on
plot(Kgrid, cons_VFI,'LineWidth',2)
%title('Consumption Policy Fucntion')
title('\sigma=5')
xlabel('k')
ylabel('g_c(k)')
legend('Linearized','VFI')
grid

% figure()
% plot(Kgrid,kprime,'LineWidth',2)
% hold on
% plot(Kgrid, kprime_VFI,'LineWidth',2)
% title('\sigma=10')
% xlabel('k')
% ylabel('g_k(k)')
% legend('Linearized','VFI')
% grid

% figure()
% plot(Kgrid,cons,'LineWidth',2)
% title('Consumption policy function')
% xlabel('k')
% ylabel('g_c(k)')
% grid
% 
% figure()
% plot(Kgrid,kprime,'LineWidth',2)
% title('Capital policy function')
% xlabel('k')
% ylabel('g_k(k)')
% grid



