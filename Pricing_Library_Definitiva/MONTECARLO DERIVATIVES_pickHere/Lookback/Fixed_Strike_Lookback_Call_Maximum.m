%% Fixed strike loockback Call on Maximum

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
DiscPayoff=exp(-r*T)*max(max(S,[],2)-K,0);
% 3. Compute the Price
[Price,~,CI]=normfit(DiscPayoff)

len1=CI(2)-CI(1)

%% 2. Antithetic Variable MC

% 1.Simulate
[S,Sav]=
% 2. Discounted Payoff
DiscPayoff=exp(-r*T)*max(max(S,[],2)-K,0);
DiscPayoffav=exp(-r*T)*max(max(Sav,[],2)-K,0);
% 3. Compute the Price
[Price,~,CI]=normfit(0.5*(DiscPayoff+DiscPayoffav))

%% 3. Control Variable MC

% 0. choose a control variable: f=ST
Ef=S0*exp( r*T ); %expected value
% 1. estimate alpha
S=assetXXXX(S0,r,par,Nsim/100,N,T);
f=S(:,end);
g=exp(-r*T)*max(max(S,[],2)-K,0);
VC=cov(f,g);
alpha=-VC(1,2)/VC(1,1);
% 2. compute the price
S=assetXXXX(S0,r,par,Nsim,N,T);
f=S(:,end);
g=exp(-r*T)*max(max(S,[],2)-K,0);
[Price,~,CI]=normfit(g+alpha*(f-Ef))