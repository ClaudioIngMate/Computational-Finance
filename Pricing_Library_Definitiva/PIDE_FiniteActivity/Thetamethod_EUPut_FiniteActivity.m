clear; close all;

% Price an EU Call Option
% Finite Activity Case (choose MERTON or KOU model)
% exploiting Theta_Method + Operator Splitting

%% Contract params, spot and risk-free rate

r=0.001;  S0=100; %market params
T=1; K=105; %contract params

%% Calibrated Parameters MERTON Model
sigma=0.6; lambda=2; mu=0.01; delta=0.2; % Merton model parameters
nu=@(y) lambda*exp(-(y-mu).^2/(2*delta^2))/sqrt(2*pi*delta^2); % Levy measure Merton 

% %% Calibrated Parameters KOU Model
% sigma=0.6; lambda=2; p=0.45; lambdaplus=50; lambdaminus=50;
% nu=@(y) p.*lambda.*lambdaplus.*exp(-lambdaplus.*y).*(y>0) ...
%          + (1-p).*lambda.*lambdaminus.*exp(-lambdaminus.*abs(y)).*(y<0);

%% Discretization params / Theta method param
N=2000; M=100; %discretization

% theta=1; % Explicit Euler (better to avoid)
% theta=0; % Implicit Euler (unconditional stable always, OK go with this)
theta=0.5; % Crank-Nickolson (unconditional stable if theta<=0.5, and payoff C^2)
% it has an improvement on the error O(dt^2+dx^2)

%% Grid
xmin=log(0.1*S0/S0);        % Smin=0.1*S0
xmax=log(4);                % Smin=4*S0
x=linspace(xmin,xmax,N+1); dx=x(2)-x(1);
dt=T/M;
%% Compute alpha and lambda (even if lambda is known, is a check)

[alpha,lambda_num,lb,ub]=levy_integral(nu,N);
Integral_Error=lambda-lambda_num % if small, good approximation

%% Matrix
% w.r.t. B&S PDE, I add lambda, alpha
A=(1-theta)*(-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));
B=-1/dt+(1-theta)*(-sigma^2/(dx^2)-(r+lambda));
C=+(1-theta)*( (r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));
At=-theta*(-(r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));
Bt=-1/dt-theta*(-sigma^2/(dx^2)-(r+lambda));
Ct=-theta*( (r-sigma^2/2-alpha)/(2*dx)+sigma^2/(2*dx^2));
M1=sparse(N+1,N+1); M2=sparse(N+1,N+1);
M1(1,1)=1; M1(end,end)=1;
for i=1:N-1
     M1(i+1,[i i+1 i+2])=[A B C];
     M2(i+1,[i i+1 i+2])=[At Bt Ct];
end

%% Backward in time solution
% At maturity
V=max(K-S0*exp(x'),0); %%%%%%%%%%%%%%%%%%% Payoff at maturity (EU Put)
BC=zeros(N+1,1); %%%%%%%%%%%%%%% Boundary condition on xMAX (EU Put)

for j=M:-1:1
    % known t_j --> unknown t_{j-1}
    BC(1)=K*exp(-r*(T-(j-1)*dt))-S0*exp(xmin); %%%%%% Boundary condition xMIN
    % Compute the Integral at time t_j
    I=levy_integral2(lb,ub,x,V,nu,S0,K,exp(-r*(T-j*dt)));
    % Solve the linear system
    V=M1\(M2*V+BC-I);
end

plot(x,V); title('Solution');
price=interp1(x,V,0,'spline') % log(S0/S0)=0

%% Auxiliary functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha,lambda_num,lb,ub]=levy_integral(nu,N)
% 1. Integral domain truncation
step=0.5; tol=10^-10;
lb=-step;
while nu(lb)>tol
    lb=lb-step;
end
ub=+step;
while nu(ub)>tol
    ub=ub+step;
end
% 2. Quadrature
N_q=2*N; y=linspace(lb,ub,N_q);
figure; plot(y,nu(y)); title('Levy measure')
alpha=trapz(y,(exp(y)-1).*nu(y));
lambda_num=trapz(y,nu(y));
end

function I=levy_integral2(lb,ub,x,V,nu,S0,K,disc)
I=zeros(size(V));
N_q=length(x);
y=linspace(lb,ub,N_q)'; dy=y(2)-y(1); nu_y=nu(y);
w=ones(N_q,1); w(1)=0.5; w(end)=0.5; % trapezoidal quadrature
% I(1)=I(end)=0
for i=2:length(I)-1
    I(i)=sum( w.*V_f(x,V,x(i)+y,S0,K,disc).*nu_y )*dy;
end
end

function v=V_f(x,V,y,S0,K,disc)
v=zeros(size(y));
%% y<=xmin     BOUNDARY CONDITION on xMIN
index=(y<=x(1));
v(index)=K*disc-S0*exp(y(index)); %%%%%%%% Boundary condition xMIN (EU Put)
%% y>= xmax    BOUNDARY CONDITION on xMAX
index=(y>=x(end));
v(index)=0; %%%%%%%%%% Boundary condition xMAX (EU Put)
%% xmin <y <xmax
index=find((y>x(1)).*(y<x(end)));
v(index)=interp1(x,V,y(index)); %linear interpolation
end
