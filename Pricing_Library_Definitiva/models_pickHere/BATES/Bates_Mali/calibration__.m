%% INITIALIZATION AND DATA
clear; close all;clc

TimeToMaturity = yearfrac('29-Oct-21','16-Sep-22',0); 
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
        56.63	95
];
[temp,index]=sort(Data(:,2));
Data=Data(index,:);
%Data=[Strike,Time to Maturity,Market Price]
Data=[Data(:,2), TimeToMaturity*ones(size(Data(:,2))), Data(:,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spot=149.80;   rf=0.13/100;
strike=Data(:,1); 
maturity=Data(:,2);
pmkt=Data(:,3); %market price
ImpVolatilityMARKET = blsimpv...
    (spot, strike, rf, maturity, pmkt);



%% Calibrating the BATES
%par=lsqnonlin(@(x) priceBATES(x,spot,strike,rf,maturity,pmkt),...
%  [0.3 0.3 0 0.3 0.3 2 0.3 0.3],[0.1 0 -1 0 0.1 1 0.1 0.1], [1 2 1 2 1 50 5 1])
par=lsqnonlin(@(x) priceBATES_IV(x,spot,strike,rf,maturity,pmkt),...
  [0.3 0.3 0 0.3 0.3 2 0.3 0.3],[0.1 0 -1 0 0.1 1 0.05 0.05], [1 2 0 2 1 50 5 1])
[Error,PriceBATES]=priceBATES(par,spot,strike,rf,maturity,pmkt);

figure
plot(PriceBATES,'sm')
hold on
plot(pmkt,'+b')
legend('Calibrated Price (BATES)','Market Price')

ImpVolatilityBATES = blsimpv(spot,...
    strike, rf, maturity, PriceBATES);

figure
plot(ImpVolatilityBATES,'sm')
hold on
plot(ImpVolatilityMARKET,'+b')
legend('Calibrated ImpVol (BATES)','Market ImpVol')

