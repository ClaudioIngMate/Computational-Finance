%% ---------------- Computational Finance Exam --------------------- %% 

% Name and Surname:     Claudio Manzoni       
% Polimi e-mail:        claudio1.manzoni@mail.polimi.it 
% Codice Persona:       10580671
% Matricola:            970219

%% Initialization and Data
clear; close all; clc

% Market parameters
S0=149.80;          % spot price
r=0.13/100;         % risk-free rate

% Options data
TimeToMaturity = 1;
% TimeToMaturity=yearfrac('29-Oct-21','16-Sep-22',0); 
%[Price,Strike]
Data=[  30.5	125
        11.4	160
        15.5	150
        5.88	180
        27.15	130
        20.75	140
        8.21	170
        3.25	200
        7.00	175
        43.18	110
        15.65	150
        31.00	125
        23.95	135
        52.03	100
        56.63	95];
    
% sort data according to the strikes (if not already sorted)
[temp,index]=sort(Data(:,2));
Data=Data(index,:);

%%% CALL option Data
%Data=[Strike,Time to Maturity,Market Price]
Data=[Data(:,2), TimeToMaturity*ones(size(Data(:,2))), Data(:,1)];

%%%%% for calibration code I need these 3 vectors: %%%%%%%%%%%%%%%%%%%%%%%%

strike=Data(:,1);           % vector of strike prices
maturity=Data(:,2);         % vector of maturities
pmkt=Data(:,3);             % vector of market prices

%%% if PUT option Data --> transform them in Call option data using
%%% Put-call parity ------  Call =  Put+ S0 - K*exp(-r*(T-t)) 
%pmkt=pmkt + S0-K*exp(-r*maturity);

ImpVolatilityMARKET = blsimpv(S0, strike, r, maturity, pmkt);

%% Calibrating the BATES

% Params to calibrate: [theta, csi, rho, eta, V0, mu, delta, lambda]

% ---------------- "Heston" parameters ---------------- %
% theta:       vol-of-var
% csi:         mean reversion speed
% rho:         correlation
% eta:         mean
% V0:          starting value of the variance (unknwown)

% ------------------ jump parameters ------------------ %
% mu:          proportional to the mean of the jump size
% delta:       volatility of the jump size
% lambda:      intensity of the Poisson process

% par=lsqnonlin(@(x) priceBATES_IV(x,S0,strike,r,maturity,pmkt),...
%   [0.3 0.3 0 0.3 0.3 2 0.3 0.3],[0.1 0 -1 0 0.1 1 0.05 0.05], [1 2 0 2 1 50 5 1])
par=lsqnonlin(@(x) priceBATES_IV(x,S0,strike,r,maturity,pmkt),...
  [0.3 0.3 0 0.3 0.3 2 0.3 0.3],[0.1 0 -1 0 0.1 1 0.05 0.05], [1 2 0 2 1 50 5 1])

[Error,PriceBATES]=priceBATES(par,S0,strike,r,maturity,pmkt);

calibration_error=sum(Error);

% Plot of Prices
figure
subplot(1,2,1);
plot(PriceBATES,'sm','Markersize', 7, 'Linewidth', 1.5)
hold on
plot(pmkt,'+b','Markersize', 7, 'Linewidth', 1.5)
legend('Calibrated Price (BATES)','Market Price')
title('Prices fitting'); grid on;

ImpVolatilityBATES = blsimpv(S0,strike, r, maturity, PriceBATES);

% Plot of Implied Volatilities
subplot(1,2,2);
plot(ImpVolatilityBATES,'sm','Markersize', 7, 'Linewidth', 1.5)
hold on
plot(ImpVolatilityMARKET,'+b','Markersize', 7, 'Linewidth', 1.5)
legend('Calibrated ImpVol (BATES)','Market ImpVol')
title('ImpliedVol fitting'); grid on;

sgtitle('First calibration');

