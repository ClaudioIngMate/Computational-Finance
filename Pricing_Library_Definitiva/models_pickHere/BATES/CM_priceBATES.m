function P_i = CM_priceBATES(param,S0,K_i,r,T)
 
% Price of an EU Call using Carr & Madan algorithm
% Bates model
% 
% INPUT         K_i       strikes (could be a vector)
% OUTPUT        P_i       prices for the K_i

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
    P_i = interp1(K,P,K_i,'spline');

    %% Plot

%     figure
%     plot(K,P, 'r') ;
%     hold on
%     axis([0  2*S0  0  S0]) ;
%     xlabel('strike') ;
%     ylabel('option price') ;
    
end

function f = CharExpBates(param,T,u)
% Characteristic Exponent  of Bates model

    % Model parameters
    theta=param(1);         % vol-of-var
    csi=param(2);           % mean reversion level
    rho=param(3);           % correlation
    eta_=param(4);          % mean
    V0=param(5);            % starting value of the variance (calibrated)
    mu=param(6);            % linked to the mean of jump size
    delta=param(7);         % volatility of the jump size
    lambda=param(8);        % intensity of the Compound Poisson process
    % Remark: log(1+jumpzise) ~ N ( log(1+mu)-0.5delta^2, delta^2 )
    
    d=sqrt((rho*theta*1i*u-csi).^2+theta^2*(1i*u+u.^2));
    g=(csi-rho*theta*u*1i-d)./(csi-rho*theta*u*1i+d);
    B=eta_*csi/theta^2*((csi-rho*theta*1i*u-d)*T-2*log((1-g.*exp(-d*T))./(1-g)));
    C=V0/theta^2*(csi-rho*theta*1i*u-d).*(1-exp(-d*T))./(1-g.*exp(-d*T));
    D=-lambda*mu*1i*u*T+lambda*T*((1+mu).^(1i*u).*exp(0.5*delta^2*1i*u.*(1i*u-1))-1);
    y = B + C + D;
    f = exp(y);
    
end
