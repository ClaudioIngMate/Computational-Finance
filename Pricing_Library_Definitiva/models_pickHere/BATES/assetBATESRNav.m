function [S, Sav] = assetBATESRNav(S0,r,par,Nsim,N,T)

% Implementation of QE scheme for Heston model
% to simulate the underlyng and its antithetical version

% "Efficient Simulation of the Heston Stochastic Volatility Model", Andersen

%% Parameters

theta=par(1);       % vol-of-var
csi=par(2);         % mean reversion speed
rho=par(3);         % correlation
eta=par(4);         % mean
V0=par(5);          % starting value of the variance (calibrated)

% jumps parameters
mu=par(6);          % linked to the mean of jump size
delta=par(7);       % volatility of the jump size
lambda=par(8);      % intensity of the Compound Poisson process
% Remark: log(1+jumpzise) ~ N ( log(1+mu)-0.5delta^2, delta^2 )

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

% % fix drift to obtain risk-neutrality
% CharExp = @(v) lambda*(exp(-delta^2*v.^2/2+1i*mu*v)-1); % jumps characteristic exponent without drift
% drift_rn = - CharExp(-1i);                                               % drift chosen under the risk neutral measure
% k0 = (r+drift_rn)*dt-rho*csi*eta*dt/theta;

k0 = (r-lambda*mu)*dt-rho*csi*eta*dt/theta; % risk neutral drift compensator for jumps
k1 = gamma1*dt*(csi*rho/theta-0.5)-rho/theta;
k2 = gamma2*dt*(csi*rho/theta-0.5)+rho/theta;
k3 = gamma1*dt*(1-rho^2);
k4 = gamma2*dt*(1-rho^2);

Z = randn(Nsim, N); % Gaussians

NT=icdf('Poisson',rand(Nsim,1),lambda*T); % number of jumps up to time T (simulation)

for j=1:Nsim
    JumpTimes=sort(T*rand(NT(j),1));
    % add the diffusion component
    for i=1:N 
    S(j,i+1)=exp(log(S(j,i))+k0+k1*V(j,i)+k2*V(j,i+1)+sqrt(k3*V(j,i)+k4*V(j,i+1)).*Z(j,i));
    Sav(j,i+1)=exp(log(Sav(j,i))+k0+k1*V(j,i)+k2*V(j,i+1)-sqrt(k3*V(j,i)+k4*V(j,i+1)).* Z(j,i));
    % add the jump part if exists a jump in ((i-1)*dt, i*dt)
        for l=1:NT(j)
            if( (JumpTimes(l)>(i-1)*dt) && ( JumpTimes(l)<=i*dt ) )
                Z_jump=randn;
                
%                 Y = mu + delta * Z_jump;
%                 Yav = mu - delta * Z_jump;
                
                Y=exp(log(1+mu)-0.5.*delta.^2 +delta.*Z_jump)-1;
                Yav=exp(log(1+mu)-0.5.*delta.^2 -delta.*Z_jump)-1;
                
                S(j,i+1)=S(j,i+1)*(1+Y);
                Sav(j,i+1)=Sav(j,i+1)*(1+Yav);
%                   S(j,i+1)=S(j,i+1)*exp(Y);
%                   Sav(j,i+1)=Sav(j,i+1)*exp(Yav);
            end
        end
    end
end

end