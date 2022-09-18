function [S,Sav]=assetVGRNav(S0,r,par,Nsim,N,T)

% Simulation of Variance Gamma underlying and its antithetical version
% S(t)=S0*exp(rt+X(t))
% X(t) is a VG in the Risk-neutral measure

% Model parameters
sigma=par(1);       % volatility of the diffusion component
theta=par(2);       % drift of the diffusion component
k=par(3);           % variance of the subordinator

%Initialization
X=zeros(Nsim,N+1); % X = rt+X(t)
Xav=X;
dt=T/N;

% Risk - neutral characteristic exponent
charexp=@(u) -1/k*log(1+u.^2*sigma.^2*k/2-1i*theta*k*u); % without drift
drift=r-charexp(-1i); % to be risk-neutral

a=dt/k;

for i=1:N
    
    %simulate the subordinator
    ds=icdf('gamma',rand(Nsim,1),a,1);
    ds=k*ds;
    
    %simulate the Wiener increment
    Z=randn(Nsim,1);
    
    % add diffusion component
    X(:,i+1)=X(:,i)+ drift*dt+theta*ds+sigma*sqrt(ds).*Z;
    Xav(:,i+1)=Xav(:,i)+ drift*dt+theta*ds-sigma*sqrt(ds).*Z;
    
end

% from logreturn to price
S=S0*exp(X); Sav=S0*exp(Xav);

end
