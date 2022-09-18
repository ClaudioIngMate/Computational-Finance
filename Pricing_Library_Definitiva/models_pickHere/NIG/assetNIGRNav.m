function [S,Sav]=assetNIGRNav(S0,r,par,Nsim,N,T)

% Simulation of NIG model for the underlying and its antithetical version
% S(t) = S0*exp(rt + X(t))
% X(t) NIG process in the risk-neutral measure

% Model parameters
sigma=par(1);       % volatility of the diffusion component
theta=par(2);       % drift of the diffusion component
k=par(3);           % variance of the subordinator

% Initialization
X=zeros(Nsim,N+1); % X = rt+X(t)
Xav=X;
dt=T/N;

% Risk - neutral characteristic exponent
charexp=@(u) 1/k*(1-sqrt(1+u.^2*sigma.^2*k-2*1i*theta*k*u)); % without drift
drift=r-charexp(-1i); % to be risk-neutral

mu=dt;
lambda=dt^2/k;

for i=1:N
    
    n=randn(Nsim,1); U=rand(Nsim,1);
    y=n.^2;
    x1=mu+mu^2*y/(2*lambda)-mu/(2*lambda)*sqrt(4*mu*lambda*y+mu^2*y.^2);
    ds=x1.*(U<=mu./(x1+mu))+mu^2./x1.*(U>mu./(x1+mu));
    
    Z=randn(Nsim,1);
    
    % add diffusion component
    X(:,i+1)=X(:,i)+ drift*dt+theta*ds+sigma*sqrt(ds).*Z;
    Xav(:,i+1)=Xav(:,i)+ drift*dt+theta*ds-sigma*sqrt(ds).*Z;
    
end

% From logreturns to spot price
S=S0*exp(X);Sav=S0*exp(Xav);

end
