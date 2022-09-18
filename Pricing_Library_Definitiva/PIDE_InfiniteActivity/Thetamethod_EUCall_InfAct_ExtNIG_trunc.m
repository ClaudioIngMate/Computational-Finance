clear; close all;

% Price an EU Call Option
% Infinite Activity case: in the extended NIG model
% Theta_Method + Operator Splitting 
% Asmussen-Rosinsky truncation: approximate a IA process with a FA process, truncating
% the jumps and approximating the truncated jumps with a Brownian Motion

%% Contract params, spot and risk-free rate

r=0.001; S0=100; % market params
T=1; K=100; %contract params

%% Discretization params / Theta method param
N=2000; M=200; %discretization
epsilon=0.2;

% theta=1; % Explicit Euler (better to avoid)
% theta=0; % Implicit Euler (unconditional stable always, OK go with this)
theta=0.5; % Crank-Nickolson (unconditional stable if theta<=0.5, and payoff C^2)
% it has an improvement on the error O(dt^2+dx^2)

%% Calibrated Parameters Extended NIG Model

sigmaGBM=0.6; thetaNIG=0.04; sigmaNIG=0.2; kNIG=0.1;
% Lèvy measure:
A = thetaNIG/sigmaNIG^2;
B = sqrt(thetaNIG^2+sigmaNIG^2./kNIG)/sigmaNIG^2;
C = sqrt(thetaNIG^2+sigmaNIG^2/kNIG)/(pi*sigmaNIG*sqrt(kNIG));
nu=@(y) C.*exp(A.*y).*besselk(1,B.*abs(y))./(abs(y));

%% Grid
xmin=log(0.1*S0/S0);            % Smin=0.1*S0
xmax=log(4);                    % Smax=4*S0
x=linspace(xmin,xmax,N+1); dx=x(2)-x(1);
dt=T/M;

%% Compute alpha and lambda

[alpha,lambda,sigma_eps,lb,ub]=levy_integral(nu,N,epsilon)
% I have to use the TRUNCATED LEVY MEASURE:
nu=@(y) (C.*exp(A.*y).*besselk(1,B.*abs(y))./(abs(y))).*(abs(y)>epsilon);
figure
y=linspace(lb,ub,N);
plot(y,nu(y));

%% Matrix
sigma=sqrt(sigmaGBM^2+sigma_eps^2); % independent brownian motions
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
V=max(S0*exp(x')-K,0); %%%%%%%%% Payoff at maturity
BC=zeros(N+1,1); %%%%%%%%%%%%% Boundary conditions
for j=M:-1:1
    % known t_j --> unknown t_{j-1}
    BC(end)=S0*exp(xmax)-K*exp(-r*(T-(j-1)*dt));
    % Compute the Integral at time t_j
    I=levy_integral2(lb,ub,x,V,nu,S0,K,exp(-r*(T-j*dt)));
    % Solve the linear system
    V=M1\(M2*V+BC-I);
end
plot(x,V); title('Solution');
price=interp1(x,V,0,'spline') % log(S0/S0)=0

%% Auxiliary functions:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha,lambda,sigma_eps,lb,ub]=levy_integral(nu,N,epsilon)
% 1. Integral domain truncation
step=0.5; tol=10^-10;
lb=-max(step,epsilon);
while nu(lb)>tol
    lb=lb-step;
end
ub=+max(step,epsilon);
while nu(ub)>tol
    ub=ub+step;
end
% 2. Quadrature
N_q=N; y=linspace(lb,-epsilon,N_q);
figure; plot(y,nu(y)); title('Levy measure')
alpha=trapz(y,(exp(y)-1).*nu(y));
lambda=trapz(y,nu(y));
N_q=N; y=linspace(epsilon,ub,N_q);
hold on; plot(y,nu(y)); 
alpha=alpha+trapz(y,(exp(y)-1).*nu(y));
lambda=lambda+trapz(y,nu(y));
y=linspace(-epsilon,epsilon,2*N); % N even
sigma_eps=sqrt(trapz(y,y.^2.*nu(y))); % definition of A-R theorem
end

function I=levy_integral2(lb,ub,x,V,nu,S0,K,disc)
I=zeros(size(V));
N_q=length(x)-1; %even number to avoid nu(0)
y=linspace(lb,ub,N_q)'; dy=y(2)-y(1); nu_y=nu(y);
w=ones(N_q,1); w(1)=0.5; w(end)=0.5; % trapezoidal quadrature
% I(1)=I(end)=0
for i=2:length(I)-1
    I(i)=sum( w.*V_f(x,V,x(i)+y,S0,K,disc).*nu_y )*dy;
end
end

function v=V_f(x,V,y,S0,K,disc)
v=zeros(size(y));
%% y<=xmin                  %%%%%% BOUNDARY CONDITION (extend...)
%index=(y<=x(1));
%v(index)=0; %%%%% since it is EU Call option
%% y>= xmax                 %%%%%% BOUNDARY CONDITION (extend...)
index=(y>=x(end));
v(index)=S0*exp(y(index))-K*disc;   %%%%% since it is EU Call option
%% xmin<y<xmax              
index=find((y>x(1)).*(y<x(end)));
v(index)=interp1(x,V,y(index));
end
