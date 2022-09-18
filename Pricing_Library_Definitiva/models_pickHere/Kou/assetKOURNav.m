function [S,Sav]=assetKOURNav(S0,r,par,Nsim,N,T)

% Simulation of KOU model for the underlying and its antithetical version
% S(t)=S0*exp(rt+X(t))
% X(t) is a Kou process in the RiskNeutral measure

% Model parameters
sigma=par(1);       % volatility of the diffusion component
lambda=par(2);      % intensity of Poisson Process
p=par(3);           % probability of positive jumps
lambdaplus=par(4);  % law of positive jumps
lambdaminus=par(5); % law of negative jumps

% Initialization
X=zeros(Nsim,N+1); %X=rt+X(t)
Xav=X;
dt=T/N;

% Number of jumps
NT=icdf('Poisson',rand(Nsim,1),lambda*T);

% Risk - neutral characteristic exponent
charexp=@(u) -sigma^2*u.^2/2+...
    1i*u*lambda.*(p./(lambdaplus-1i*u)-(1-p)./(lambdaminus+1i*u)); % without drift
drift=r-charexp(-1i); % to be risk-neutral

for j=1:Nsim
    JumpTimes=sort(T*rand(NT(j),1));
    Z=randn(N,1);
    for i=1:N
        
      % add diffusion component
      X(j,i+1)=X(j,i)+drift*dt+sigma*sqrt(dt)*Z(i);
      Xav(j,i+1)=Xav(j,i)+drift*dt-sigma*sqrt(dt)*Z(i);
      
      % add jump part -> only if a jump is present in ( (i-1)dt,idt ]
      for l=1:NT(j)
          if JumpTimes(l)>(i-1)*dt && JumpTimes(l)<=i*dt
              sim_p=rand;
              if sim_p<p %positive jump
                    Y=icdf('exp',rand,1/lambdaplus);
              else %negative jump
                    Y=-icdf('exp',rand,1/lambdaminus);
              end
              X(j,i+1)=X(j,i+1)+Y;
              Xav(j,i+1)=Xav(j,i+1)+Y;
          end
      end
    end
end

% From logreturn to underlying spot price
S=S0*exp(X);
Sav=S0*exp(Xav);

end