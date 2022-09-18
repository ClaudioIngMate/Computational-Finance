function [P_i] = CM_priceHESTON(x,S0,K_i,r,T)

% Price of a European Call exploiting Carr & Madan algorithm
% Heston model

% INPUT         K_i     strikes (could be a vector)
% OUTPUT        P_i     prices for the K_i

%% Discretization Parameters

Npow=15; N=2^Npow;  % discretization points
A=600;              % upper bound

%% FFT Grid (for dv and dk)

% dv
eta=A/N; v = eta*(0:N-1);   % integral domain grid (only > 0)
v(1)=1e-22;                 % correction term: could not be equal to zero (otherwise NaN)                    

% dk
lambda=2*pi/(N*eta); k=-lambda*N/2+lambda*(0:N-1);  % log-strike grid  
K = S0 * exp(k); % strike

%% FT{Z(k)} using Carr&Madan formula

tr_fou = trasf_fourier(r,x,T,v);                % Fourier transform of z_k

%% Prices computing the integral using trapezoidal rule

w = [0.5  ones(1,N-2)  0.5];                    % Trapezoidal rule 
h = exp(1i*(0:N-1)*pi).*tr_fou.*w*eta;
P = S0 * real( fft(h)/pi + max(1-exp(k-r*T),0));% prices

%%  Output: filter prices grid + interpolation

index=find( (K>0.1*S0 & K<3*S0) );      % delete nonsense strikes
K=K(index); P=P(index);
P_i = interp1(K,P,K_i,'spline');  % option price interpolating on the pricing grid

%%%% -------------------- PLOT -------------------- %%%% 
% PLOT
% figure
% plot(K,P, 'r') ;
% hold on
% axis([0  2*S0  0  S0]) ;
% xlabel('strike') ;
% ylabel('option price') ;
%%%% ----------------------------------------------- %%%%

%% Auxiliary functions

function fii = trasf_fourier(r,x,T,v)
fii = (exp(1i*r*v*T)).*( (characteristic_func(x,T,v-1i)-1)./(1i*v.*(1+1i*v))  );
end

function f = characteristic_func(x,T,u)

% Model parameters
theta=x(1);         % vol-of-var
csi=x(2);           % mean reversion speed
rho=x(3);           % correlation
eta_=x(4);          % mean 
V0=x(5);            % starting point of the variance
    
alfa = -.5*(u.*u + u*1i);
beta = csi - rho*theta*u*1i;
epsilon2 = theta * theta;
gamma = .5 * epsilon2;

D = sqrt(beta .* beta - 4.0 * alfa .* gamma);

bD = beta - D;
eDt = exp(- D * T);

G = bD ./ (beta + D);
B = (bD ./ epsilon2) .* ((1.0 - eDt) ./ (1.0 - G .* eDt));
psi = (G .* eDt - 1.0) ./(G - 1.0);
A_ = ((csi * eta_) / (epsilon2)) * (bD * T - 2.0 * log(psi));

y = A_ + B*V0;

f = exp(y);

end

end