function [Price] = CM_priceBATES(param,S0,Strike,r,T)
%
% Price of an Eu Call using Carr-Madan algorithm under Bates model.
%
% INPUT: 
%  - param: Bates model parameters:
%     epsilon: vol-of-vol  
%     kappa: mean reversion speed    
%     rho: correlation
%     theta: mean reversion level  
%     V0: initial variance
%     lambda: jumps intensity
%     kbar, delta: log(1+jumpzise) ~ N ( log(1+kbar)-0.5delta^2, delta^2 )
%  - S0: spot price
%  - Strike: strike (could be a vector)
%  - r: risk-free interest rate
%  - T: maturity of the option
% 
% OUTPUT: 
% P_i = prices (could be a vector)
%

    %% Grid

    Npow = 20; N = 2^(Npow); % grid point
    A = 600; % upper bound
    eta = A/N; lambda = 2*pi/(N*eta); 
    k = -lambda * N/2  + lambda *(0:N-1); % log-strike grid  
    K = S0 * exp(k); % strike 
    v = eta*(0:N-1); v(1)=1e-22; %correction term: could not be equal to zero (otherwise NaN)
    
    %% Pricing

    % Fourier transform of z_k
    CharFun = (exp(1i*r*v*T)).*( (CharExpBates(param,T,v-1i)-1)./(1i*v.*(1+1i*v))  );

    % Trapezoidal rule
    w = [0.5  ones(1,N-2)  0.5];  
    h = exp(1i*(0:N-1)*pi).*CharFun.*w*eta;
    P = S0 * real( fft(h)/pi + max(1-exp(k-r*T),0)); % prices
    
    % delete too small and too big strikes
    index=find( (K>0.1*S0 & K<3*S0) );
    K=K(index); P=P(index);

    % Interpolation
    Price = interp1(K,P,Strike,'spline');

    %% Plot
    % Da sopprimere nella calibration!!!

%     figure
%     plot(K,P, 'r') ;
%     hold on
%     axis([0  2*S0  0  S0]) ;
%     xlabel('strike') ;
%     ylabel('option price') ;
    
   
end

function f = CharExpBates(param,T,u)

    epsilon=param(1); kappa=param(2); rho=param(3); theta=param(4);
    V0=param(5); lambda=param(6); kbar=param(7); delta=param(8); 
    
    d=sqrt((rho*epsilon*1i*u-kappa).^2+epsilon^2*(1i*u+u.^2));
    g=(kappa-rho*epsilon*u*1i-d)./(kappa-rho*epsilon*u*1i+d);
    
    B=theta*kappa/epsilon^2*((kappa-rho*epsilon*1i*u-d)*T-2*log((1-g.*exp(-d*T))./(1-g)));
    C=V0/epsilon^2*(kappa-rho*epsilon*1i*u-d).*(1-exp(-d*T))./(1-g.*exp(-d*T));
    D=-lambda*kbar*1i*u*T+lambda*T*((1+kbar).^(1i*u).*exp(0.5*delta^2*1i*u.*(1i*u-1))-1);
    
    y = B + C + D;
    
    f = exp(y);
end
