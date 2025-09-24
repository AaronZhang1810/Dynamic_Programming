%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firm price adjustment problem
% Wrapper for model solution and simulation
% Firm_Pricing_Wrapper.m
%
% S. Terry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

disp('%%%%%%%%%%%%% Solving the firm price adjustment problem')
disp(' ')

%%%%%%%%%% set model parameters

epsilon = 2; %demand elasticity
c = 0.75; %quadratic price adjustment cost
beta = 0.96; %discount rate
rho = 0.85; %persistence of firm marginal cost m
sigma = 0.05; %standard deviation of marginal cost innovation

%%%%%%%%%% set solution parameters

%grid parameters
pnum = 1000; %number of price grid points
pmin = 1; %lowest grid point for price
pmax = 3.5; %highest grid point for price
mnum = 35; %number of marginal cost shocks grid points
mnstdev = 3; %number of standard deviations to cover around steady state marginal cost
statenum = mnum*pnum; %total number of state grid points

%solution parameters
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

%some formatting stuff
lwidnum = 2;
fsizenum = 12;

%%%%%%%%%% set up the grids, discretization of marg cost, payoff matrices
Firm_Pricing_Setup;

%%%%%%%%%% implement simple VFI
%Firm_Pricing_VFI;

%%%%%%%%%% implement VFI with Howard Improvement
Firm_Pricing_VFI_Howard;

%%%%%%%%%% implement VFI with MacQueen-Porteus Bounds
Firm_Pricing_VFI_MQP;

%%%%%%%%%% simulate the model
Firm_Pricing_Simulate;

%%%%%%%%%% simulate the model with nonstochastic simulation
Firm_Pricing_Simulate_Dist;
