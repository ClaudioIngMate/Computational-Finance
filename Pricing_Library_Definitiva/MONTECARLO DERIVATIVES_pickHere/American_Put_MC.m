%% American Put option with MonteCarlo (L&S algorithm)

% Market parameters
S0=;        % underlying spot price
r=;         % risk-free interest rate

% Contract Parameters
K=;
T=;

% Model parameters vector -----> same order of ASSET_SIMULATION code!!
par=[ ,  ]

Nsim=1e5;% Number of simulations

N=round(T*12); % monthly monitoring
N=round(T*52); % weekly monitoring
N=round(T*252); % business day monitoring
N=round(T*365); % daily monitoring

dt=T/N;

%% Pricing

% Asset simulation
SPaths=assetXXXXXX(....

% INITIALIZATION
ExerciseTime=M*ones(Nsim,1); %initialize at M
CashFlows=max(0,K-SPaths(:,M)); %payoff


% BACKWARD IN TIME PROCEDURE
% (we consider only the InMoney case, thus option will be exercised)
for step=M-1:-1:1
    
    InMoney=find(SPaths(:,step)<K);
    S=SPaths(InMoney,step);
    
    %%%%---------------- Regression ----------------%%%%
    % basis functions = [1,S,S^2]
    RegrMat=[ones(length(S),1), S, S.^2];
    YData=CashFlows(InMoney).*exp(-r*dt*(ExerciseTime(InMoney)-step));
    alpha=RegrMat\YData;
    
    %%%%----------------- IV and CV ----------------%%%%
    IV=K-S; %intrinsic value
    CV=RegrMat*alpha; %continuation value (its approximation)
    
    %%%%-------------- Early exercise --------------%%%%
    % Paths with early exercise at time step
    Index=find(IV>CV);
    ExercisePaths=InMoney(Index); %index of simulation where it is optimal to early exercise at time step
    
    % Update Cashflows
    CashFlows(ExercisePaths)=IV(Index);
    % Update Exercise Time
    ExerciseTime(ExercisePaths)=step;
end

[Price,~,CI]=normfit(CashFlows.*exp(-r*dt*ExerciseTime))

