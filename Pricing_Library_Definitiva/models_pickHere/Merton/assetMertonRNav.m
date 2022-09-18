function [S,Sav]=assetMertonRNav(S0,r,par,Nsim,N,T)

% Simulation of MERTON model for the underlying and its antithetical version
% S(t) = S0* exp (rt + X(t))
% X(t) is a Merton process in the risk-neutral measure

% Model parameters
sigma=par(1);       % volatility of the diffusion component
lambda=par(2);      % intensity of Poisson process
mu=par(3);          % mean value of the jump size
delta=par(4);       % standard deviation of the jump size

% Initialization
X=zeros(Nsim,N+1); %X=rt+X(t)
Xav=X;
dt=T/N;

% Number of jumps in [0, T]
NT=icdf('Poisson',rand(Nsim,1),lambda*T);

% Risk - neutral characteristic exponent
charexp=@(u) -sigma^2*u.^2/2+ lambda*(exp(-delta^2*u.^2/2+1i*mu*u)-1); % without drift
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
            if ( (JumpTimes(l)>(i-1)*dt) && (JumpTimes(l)<=i*dt) )
                g = randn;
                Y=mu+delta*g;
                Yav=mu-delta*g;
                X(j,i+1)=X(j,i+1)+Y;
                Xav(j,i+1)=Xav(j,i+1)+Yav;
            end
        end
    end
end

% From log-returns to underlying spot price
S=S0*exp(X);
Sav=S0*exp(Xav);

end