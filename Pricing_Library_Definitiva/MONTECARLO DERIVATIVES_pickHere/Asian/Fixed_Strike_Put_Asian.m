%% Fixed Strike Asian Put option with MonteCarlo

% Market parameters
S0=;
r=;

% Contract Parameters
K=;
T=;

% Model parameters vector -----> same order of ASSET_SIMULATION code!!
par=[ ,  ]

Nsim=1e6;% Number of simulations

N=round(T*12); % monthly monitoring
N=round(T*52); % weekly monitoring
N=round(T*252); % business day monitoring
N=round(T*365); % daily monitoring

dt=T/N;

%% 1. Classical MC

% 1. Simulate
S=
% 2. Discounted Payoff
DiscPayoff=exp(-rf*T)*max(K-mean(S,2),0);
% 3. Compute the Price
[Price,~,CI]=normfit(DiscPayoff)

len1=CI(2)-CI(1)

%% 2. Antithetic Variable MC

% 1.Simulate
[S,Sav]=
% 2. Discounted Payoff
DiscPayoff=exp(-rf*T)*max(K-mean(S,2),0);
DiscPayoffav=exp(-rf*T)*max(K-mean(Sav,2),0);
% 3. Compute the Price
[Price,~,CI]=normfit(0.5*(DiscPayoff+DiscPayoffav))

len2=CI(2)-CI(1)

%% 3. Control Variable MC

% 0. choose a control variable: f=K-mean(S,2)
dt=T/N;
%expected value of f
Ef=K-spot/(N+1)*sum( exp(rf*(0:N)*dt) );
% 1. estimate alpha
S=assetXXXX(spot,rf,par,Nsim/100,N,T);
f=K-mean(S,2);
g=exp(-rf*T)*max(f,0);
VC=cov(f,g);
alpha=-VC(1,2)/VC(1,1);
% 2. compute the price
S=assetXXXX(spot,rf,par,Nsim,N,T);
f=K-mean(S,2);
g=exp(-rf*T)*max(f,0);
[Price,~,CI]=normfit(g+alpha*(f-Ef))

len3=CI(2)-CI(1)

