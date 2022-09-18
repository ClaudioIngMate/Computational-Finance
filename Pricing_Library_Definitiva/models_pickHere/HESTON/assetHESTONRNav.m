function [S, Sav] = assetHESTONRNav(S0,r,par,Nsim,N,T)

% Implementation of QE scheme for Heston model
% to simulate the underlyng and its antithetical version

% "Efficient Simulation of the Heston Stochastic Volatility Model", Andersen

%% Parameters

theta=par(1);       % vol-of-var
csi=par(2);         % mean reversion speed
rho=par(3);         % correlation
eta=par(4);         % mean
V0=par(5);          % starting value of the variance (calibrated)

dt=T/N; 
V=zeros(Nsim,N+1);
S=zeros(Nsim,N+1);
Sav=zeros(Nsim,N+1);

% Random changes for volatility are sampled from one of two distributions,
% depending on the ratio Psi = s^2/m^2, where m & s are mean and variance
% of next volatility value, conditioned to current one.
%  Scheme 1 (Exponential)        is selected when Psi > Psi_cutoff
%  Scheme 2 (Quadratic)          is selected when 0 < Psi < Psi_cutoff

Psi_cutoff = 1.5;       % as suggested in article

%% Discetize V
V(:,1)=V0;
for i = 1:N
    
    % STEP 1 & 2  --> compute m, s, Psi
    m= eta+(V(:,i)-eta)*exp(-csi*dt);
    m2= m.^2;
    s2= V(:,i)*theta^2*exp(-csi*dt)*(1-exp(-csi*dt))/csi +eta*theta^2*(1-exp(-csi*dt))^2/(2*csi);
    s = sqrt(s2);
    Psi = (s2)./(m2);
    
    % STEP 3,4,5 --> depending on Psi, use E or Q scheme to calculate next V
    % 1. Exponential approximation for Psi > Psi_cutoff
    % The PDF of V(t+dt) is p * delta(0) + (1 - p) * (1 - exp( - beta x )
    % thus a probability mass in 0 and an exponential tail after that
    index = find( Psi > Psi_cutoff );% simulations where Psi > Psi_cutoff
    p_exp =(Psi(index)-1)./(Psi(index)+1);	% probability mass in 0 (Eq. 29)
    beta_exp = (1-p_exp)./m(index);		    % exponent of exponential density tail (Eq. 30)
    % gets x from inverse CDF applied to uniform U
    U = rand(size(index));
    V(index,i+1)=(log((1-p_exp)./(1-U))./beta_exp).*(U>p_exp );
    % 2. Quadratic approximation for  0 < Psi < Psi_cutoff
    % V(t+dt) = a( b + Zv )^2, Zv ~ N(0, 1)
    index=find(Psi<=Psi_cutoff);  % simulations where 0 < Psi < Psi_cutoff
    invPsi =1./Psi(index);
    b2_quad=2*invPsi-1+sqrt(2*invPsi).*sqrt(2*invPsi-1);  %  (Eq 27)
    a_quad=m(index)./(1+b2_quad);	              %  (Eq. 28)
    V(index,i+1)=a_quad.*(sqrt(b2_quad)+randn(size(index))).^2; % (Eq.23)
    
end

%% Discetize S

S(:, 1)   = S0;
Sav(:, 1) = S0;
gamma1 = 0.5; % => central discretization scheme
gamma2 = 0.5; % => central discretization scheme
k0 = r*dt-rho*csi*eta*dt/theta;
k1 = gamma1*dt*(csi*rho/theta-0.5)-rho/theta;
k2 = gamma2*dt*(csi*rho/theta-0.5)+rho/theta;
k3 = gamma1*dt*(1-rho^2);
k4 = gamma2*dt*(1-rho^2);

Z = randn(Nsim, N); % Gaussians
for i=1:N
    S(:,i+1)=exp(log(S(:,i))+k0+k1*V(:,i)+k2*V(:,i+1)+sqrt(k3*V(:,i)+k4*V(:,i+1)).*Z(:,i));
    Sav(:,i+1)=exp(log(Sav(:,i))+k0+k1*V(:,i)+k2*V(:,i+1)-sqrt(k3*V(:,i)+k4*V(:,i+1)).* Z(:,i));
end

end