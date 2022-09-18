function [S,Sav]=assetExtNIGRNav(S0,r,par,Nsim,N,T)

% Simulation of Extended NIG model for underlying and its antithetical version
% S(t) = exp( rt + X(t))
% X(t) is ExtNIG process under the risk-neutral measure

% Model parameters
sigma=par(1);       % volatility of the diffusion component
theta=par(2);       % drift of the diffusion component
k=par(3);           % variance of the subordinator
sigmaGBM=par(4);    % volatility of the pure Brownian motion (no time change)

% Initialization
X=zeros(Nsim,N+1); % X  = rt+X(t)
Xav=X;
dt=T/N;

% Risk-neutral characteristic exponent
charexp=@(u) -(sigmaGBM^2/2)*u.^2+...
    1/k*(1-sqrt(1+u.^2*sigma.^2*k-2*1i*theta*k*u)); % without drift
drift=r-charexp(-1i); % to be risk-neutral

mu=dt;
lambda=dt^2/k;

for i=1:N
    
    n=randn(Nsim,1); U=rand(Nsim,1);
    y=n.^2;
    x1=mu+mu^2*y/(2*lambda)-mu/(2*lambda)*sqrt(4*mu*lambda*y+mu^2*y.^2);
    ds=x1.*(U<=mu./(x1+mu))+mu^2./x1.*(U>mu./(x1+mu));
    
    Z1=randn(Nsim,1);
    Z2=randn(Nsim,1);
    
    % add diffusion component
    X(:,i+1)=X(:,i)+ drift*dt+sigmaGBM*sqrt(dt)*Z1+theta*ds+sigma*sqrt(ds).*Z2;
    Xav(:,i+1)=Xav(:,i)+ drift*dt-sigmaGBM*sqrt(dt)*Z1+theta*ds-sigma*sqrt(ds).*Z2;
    
end
    
% From logreturn to underlying spot price 
S=S0*exp(X);
Sav=S0*exp(Xav);
