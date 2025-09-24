%%%%%% Firm's Price Adjustment Problem %%%%%%
% Date: 07/2023
% Author: Jihuan Zhang

%%%%%% Setup %%%%%%
% V(p_{-1},m) = max_{p}{ PI(p,m) - AC(p,p_{-1}) + beta*E[V(p,m')|m] }
% state space: (p_{-1},m)
% p_{-1}: price in the last period
% m: marginal cost shocks: ln(m') = rho*ln(m) + sigma*u; u~N(0,1)
% PI: firm's profits: PI = (p-m)p^{-epsilon}
% AC: price adjustment costs: AC = (c/2)(p-p_{-1})^2
% beta: discount factor

%% Parameters and Initialization
clear; close all; clc;

% model parameters
epsilon = 2; % demand elasticity
c = 0.75;
beta = 0.96; % discount factor
rho = 0.85; % persistency
sigma = 0.05; % standard deviation of the marginal cost innovation
sd = 3; % number of sd in tauchen to cover around steady state marginal cost

% grid setup
pnum = 1000; % number of price grid points
pmax = 3.5; % constraint of price is [1,3.5]
pmin = 1;
% make the grids sparse at large values and dense at low values
pgrid = linspace(log(pmax),log(pmin),pnum);
pgrid = exp(pgrid)';

mnum = 35; % number of marginal cost shocks grid points
statenum = mnum*pnum; % total number of state grid points

% solution parameters
soltol = 1e-7; %tolerance on model solution
maxit = 1000; %max iterations on model solution
howardnum = 10; %number of push forward periods in Howard acceleration

%simulation parameters
firmnum = 5000; %number of firms in simulation
Terg = 50; %number of burn-in periods per firm
Tsim = 15; %number of periods for useful simulation
Ttot = Terg + Tsim+1; %total number of periods in simulation
rng("default"); %set random seed for simulation draws
minit = floor(mnum/2); %initial point for m in simulation
pinit = floor(pnum/2); %inital point for p_{-1} in simulation