% Price of a Knock & Out Call with Conv. Method
% Can be adapted for whatever Out Option
clear; close all; clc;

%% Contract and market Parameters

r=0.02;                 % risk-free rate
T=1;                    % Time to maturity
S0=100;                 % underlying spot price
K=95;                   % strike price
Ndates=round(T*12);     % monthly monitoring
Barrier_L=0.6*S0;       % Lower barrier 
Barrier_U=1.4*S0;       % Upper barrier

%% Model parameters

params=[0.3 0.3 0.3];          % NIG ----------- sigma, theta, k

% params=[0.3 0.3 0.3 0.3];     % Extended NIG -- sigma, theta, k, sigmaGBM
% params=[0.3 0.3 0.3];         % VG ------------ sigma, theta, k
% params=[0.3 0.3 0.3 0.3];     % Extended VG --- sigma, theta, k, sigmaGBM
% params=[0.3 30 0.5 4 4];      % Kou ----------- sigma, lambda, p , lambda+, lambda-
% params=[0.3 50 0 0.3];        % Merton -------- sigma, lambda, mu, delta
% params=[0.3 0.5 0 -0.1 0.3];  % Heston -------- theta, csi, rho, eta, V0

%% Conv. method to price this option

N=2^12; % number of points

[S,v] = CONV(S0, K, r, T, Ndates, N, Barrier_L, Barrier_U, params);
price=interp1(S,v,S0,'spline')
