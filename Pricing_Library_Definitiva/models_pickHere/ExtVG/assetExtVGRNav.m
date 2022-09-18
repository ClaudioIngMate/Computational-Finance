function [S,Sav]=assetExtVGRNav(S0,r,par,Nsim,N,T)

% Simulation of Extended VG model for underlying and its antithetical version
% S(t) = exp (rt + X(t))
% X(t) is an ExtVG model in the risk-neutral measure

sigma=par(1);       % volatility of the diffusion component
theta=par(2);       % drift of the diffusion component
k=par(3);           % variance of the subordinator
sigmaGBM=par(4);    % volatility of the pure Brownian motion (no time change)

% Initialization
X=zeros(Nsim,N+1); % X=rt+X(t)
Xav=X;
dt=T/N;

% Risk-neutral characteristic exponent
charexp=@(u) -(sigmaGBM^2/2)*u.^2 -1/k*log(1+u.^2*sigma.^2*k/2-1i*theta*k*u); % without drift
drift=r-charexp(-1i); % to be risk-neutral

a=dt/k;

for i=1:N
    
    %simulate the subordinator
    ds=icdf('gamma',rand(Nsim,1),a,1);
    ds=k*ds;
    
    %simulate the Wiener processes
    Z1=randn(Nsim,1);
    Z2=randn(Nsim,1);
    
    % add diffusion component
    X(:,i+1)=X(:,i)+ drift*dt+sigmaGBM*sqrt(dt)*Z1+theta*ds+sigma*sqrt(ds).*Z2;
    Xav(:,i+1)=Xav(:,i)+ drift*dt-sigmaGBM*sqrt(dt)*Z1+theta*ds-sigma*sqrt(ds).*Z2;
    
end
    
% From logreturn to underlying spot price
S=S0*exp(X);
Sav=S0*exp(Xav);

end