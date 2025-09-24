%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm price adjustment problem
% Setting up grids and discretization and payoff matrices
% Firm_Pricing_Setup.m
%
% S. Terry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%% Set up grids and discretization and payoffs')
disp(' ')
tic

%%%%%%%%%% set up loglinear grid for endogenous prices

pgrid = linspace(log(pmin),log(pmax),pnum);
pgrid = exp(pgrid)'; %convert to levels, pnum x 1

%%%%%%%%%% set up loglinear grid & discretization for exogenous marginal cost process

[mgrid, pr_mat_m] = tauchen(sigma, rho, mnstdev, mnum);


%%%%%%%%%% set up joint indexes of the state space

grid_val = zeros(statenum,2); %mnum*pnum x 2, with (i,1) = value of p, (i,2) = value of m
grid_ind = zeros(statenum,2); %mnum*pnum x 2, with (i,1) = index of p, (i,2) = index of m

%insert values
grid_val(:,1) = kron(ones(mnum,1),pgrid);
grid_val(:,2) = kron(mgrid,ones(pnum,1));

%insert indexes
grid_ind(:,1) = kron(ones(mnum,1),(1:pnum)');
grid_ind(:,2) = kron((1:mnum)',ones(pnum,1));

%%%%%%%%%% set up payoff matrices

%Rmat is mnum*pnum x pnum, with (i,j) static payoff to state i, policy j
Rmat = zeros(statenum,pnum);

%ACmat is mnum*pnum x pnum, with (i,j) price adjustment costs with state i, policy j
%repmat(x,n,m) means repeat x n times and put them into different rows, and then
%repeat the result m times and put them into different columns
ACmat = (c/2)*(repmat(pgrid',statenum,1) - repmat(grid_val(:,1),1,pnum)).^2;

%PImat is mnum*pnum x pnum, with (i,j) flow profits with state i, policy j
PImat = (repmat(pgrid',statenum,1) - repmat(grid_val(:,2),1,pnum)).*(repmat(pgrid',statenum,1).^-epsilon);

%actually compute total static payoff or reward to state i, policy j in (i,j)
Rmat = PImat - ACmat;
toc
disp(' ')
