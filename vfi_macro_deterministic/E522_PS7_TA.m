% E522 PS7: Value function iteration pratice
% Seokil Kang
%---------------
% house keeping
%---------------
clear;
%close all;
clc;

%-------------
% model setup
%-------------

% VFI setup
N = 1e3;                          % grid size choice
tol = 1e-10;                       % tolerance choice

% calibration
beta = .96;                         % discount rate
alpha = .33;                        % production share of capital

% capital grid setup
k_1 = 1e-4;                         % minimum capital
k_N = 1;                            % maximum capital

% construct uniform grid
k_grid = linspace(k_1,k_N,N)';

% value function and policy function setup
v_old = zeros(N,1);                 % initial value function guess
v_new = zeros(N,1);                 % initial value function guess
k_prime = zeros(N,1);               % policy function basket

%--------------------------
% value function iteration
%--------------------------

% initial difference
d = tol+1;      % computer requires initial d greater than tol
iter = 0;       % iteration counting


% value function iteration part
while d > tol
    
    % count up each iteration
    iter = iter + 1;
    
    % solve the value function(max part) for each grid point
    for i = 1:N
        
        % select grid point
        k = k_grid(i);
        
        % find the maximum of value function at grid i
        % if consumption choice(k^alpha-k') is negative, punish that choice with very low value(1e-99)
        val = log(max(k^alpha-k_grid,1e-99)) + beta*v_old;
        
        % record the maximum value and its index
        [v_sol, indx] = max(val);
        
        % k prime is chosen within the grid
        k_prime(i) = k_grid(indx);
        
        % value function after the contraction
        v_new(i) = v_sol;
        
    end
    
    % measure the distance between the old & new value function
    % if d is shorter than the desired tolerance, while loop will end at this point
    d = max(abs((v_new-v_old)./v_old));
    
    % if d is still larger than the tolerance level update the old value function with the new one
    v_old = v_new;
    
    % report the iteration for every 10th iteration
    if mod(iter,10) == 0        
        fprintf('current iteration = %.0f with convergence = %.5f\n', iter,d)
    end
    
end

% report that the iteration has finished
fprintf('VFI convereged at %.0f iterations\n',iter)

%------------------
% reporting result
%------------------

% guess & verify method for true value function
A = alpha/(1-alpha*beta);
B = (log(1-alpha*beta) + alpha*beta/(1-alpha*beta)*log(alpha*beta))/(1-beta);
v_true = A*log(k_grid) + B;

% true policy function for capital
k_true = alpha*beta*k_grid.^alpha;

% linear approximation version
delta = 1;
sigma = 1;

% steady state values
kstar = ((1/beta - 1 + delta)/alpha)^(1/(alpha-1));
cstar = kstar^alpha - delta * kstar;
wstar = kstar^alpha - kstar * alpha*kstar^(alpha-1);

% second order derivative for production function w.r.t k
f_double_prime_k = alpha * (alpha - 1) * kstar^(alpha-2);

% matrix system following lecture note(The Neoclassical Growth Model, p12)
A1 = [ 1 -beta*f_double_prime_k*cstar*sigma; 0 1];
A0 = [1 0; -1 1/beta];
b = [-beta*f_double_prime_k*cstar*sigma*kstar; wstar];
A = inv(A1)*A0;
Ab = inv(A1)*b;

% eigenvalue decomposition
[P, L] = eig(A);

% inverse of eigenvector matrix
P_inv = inv(P);

% constant vector of linear system
G = P_inv * Ab;

% constant term of solution wrt unstable eigenvalue (consumption)
Gc = 1/P_inv(1,1)*inv(eye(1)-L(1))*G(1);

% linear approximation of consumption policy function
clin = -P_inv(1,2)/P_inv(1,1)*k_grid + Gc;

% linear approximation of capital policy function
kplin = L(end)*k_grid + inv(P_inv(2,2)-P_inv(2,1)*P_inv(1,2)/P_inv(1,1)) * (G(2) + (L(end)-1)*P_inv(2,1)*Gc);


%% Report results

% plot settings
lw1 = 6;
lw2 = 3;
ms = 6;
mycol = {'#00539a' '#ff832b' '#198038' '#a56eff'};
ftsz = 36;

figure('name','value function iteration','color','w','WindowState','maximized')
nexttile
hold on
p1 = plot(k_grid, v_true,'color',mycol{1},'linewidth',lw1);
p2 = plot(k_grid,v_new,'o','color',mycol{2},'linewidth',lw2,'markersize',ms);
hold off
grid on
legend([p1 p2],'true value function','VFI result','location','northwest')
legend boxoff
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
xlabel('capital grid')
title('value function','FontName','Consolas','fontsize',ftsz,'FontWeight','bold');

nexttile
hold on
p1 = plot(k_grid, k_true,'color',mycol{1},'linewidth',lw1);
p2 = plot(k_grid,k_prime,'o','color',mycol{2},'linewidth',lw2,'markersize',ms);
p3 = plot(k_grid, kplin,':','color',mycol{3},'linewidth',lw2);
hold off
grid on
legend([p1 p2 p3],'true policy function','VFI result','linear approximation','location','northwest')
ticklabel = get(gca,'TickLabel');
set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
xlabel('capital grid')
title('policy function','FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
legend boxoff