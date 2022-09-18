clc, clear; close all;

% Price an EU Call Option
% Infinite Activity Case (Here: Extended VG model)
% Theta_Method + Operator Splitting

%% Contract params, spot and risk-free rate

r=0.001; S0=100; % market params
T=1; K=105; %contract params

%% Discretization params / Theta method param
N=2000; M=200; %discretization

% theta=1; % Explicit Euler (better to avoid)
%theta=0; % Implicit Euler (unconditional stable always, OK go with this)
theta=0.5; % Crank-Nickolson (unconditional stable if theta<=0.5, and payoff C^2)
% it has an improvement on the error O(dt^2+dx^2)

%% Calibrated Parameters Extended VG Model

sigmaGBM=0.6; thetaVG=0.04; sigmaVG=0.2; kVG=0.1;
% Lèvy measure:
A=thetaVG/sigmaVG^2; B=sqrt(thetaVG^2+2*sigmaVG^2/kVG)/sigmaVG^2;
nu=@(y) exp(A*y-B*abs(y))./(kVG*abs(y));
[lb,ub]=levy_integral(nu,N); % only plots Lèvy measure. IA -> NO alpha, lambda

%% Grid
xmin=log(0.1*S0/S0);            % Smin=0.1*S0
xmax=log(4);                    % Smax=4*S0
x=linspace(xmin,xmax,N+1); dx=x(2)-x(1);
dt=T/M;

%% Matrix
% same as B&S PDE, since we do not split the integral
A=(1-theta)*(-(r-sigmaGBM^2/2)/(2*dx)+sigmaGBM^2/(2*dx^2));
B=-1/dt+(1-theta)*(-sigmaGBM^2/(dx^2)-r);
C=+(1-theta)*( (r-sigmaGBM^2/2)/(2*dx)+sigmaGBM^2/(2*dx^2));
At=-theta*(-(r-sigmaGBM^2/2)/(2*dx)+sigmaGBM^2/(2*dx^2));
Bt=-1/dt-theta*(-sigmaGBM^2/(dx^2)-r);
Ct=-theta*( (r-sigmaGBM^2/2)/(2*dx)+sigmaGBM^2/(2*dx^2));
M1=sparse(N+1,N+1); M2=sparse(N+1,N+1);
M1(1,1)=1; M1(end,end)=1;
for i=1:N-1
     M1(i+1,[i i+1 i+2])=[A B C];
     M2(i+1,[i i+1 i+2])=[At Bt Ct];
end

%% Backward in time solution
% At maturity
V=max(S0*exp(x')-K,0); %%%%%%%%%% Payoff at maturity
BC=zeros(N+1,1); %%%%%%%%%%%% Boundary condition

for j=M:-1:1
    % known t_j --> unknown t_{j-1}
    BC(end)=S0*exp(xmax)-K*exp(-r*(T-(j-1)*dt)); %%%%%%%%% Boundary condition
    % Compute the Integral at time t_j
    I=levy_integral2(lb,ub,x,V,nu,S0,K,exp(-r*(T-j*dt)));
    % Solve the linear system
    V=M1\(M2*V+BC-I);
end

%% Output

plot(x,V); title('Solution');
price=interp1(x,V,0,'spline') % log(S0/S0)=0

%% Auxiliary functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lb,ub]=levy_integral(nu,N)
% 1. Integral domain truncation
step=0.5; tol=10^-10;
lb=-step;
% find the lower bound as the value s.t. the Lèvy measure is > tolerance
while nu(lb)>tol
    lb=lb-step;
end
ub=+step;
% find the upper bound as the value s.t. the Lèvy measure is > tolerance
while nu(ub)>tol
    ub=ub+step;
end
% 2. Plot
y=linspace(lb,ub,N);% domain
figure; plot(y,nu(y)); title('Levy measure')
end

function I=levy_integral2(lb,ub,x,V,nu,S0,K,disc)
% Computes the "unknown" integral (nonlinear part)

I=zeros(size(V)); % I(1)=I(end)=0
N_q=length(x)-1;% Number of quadrature points
% Remark: avoid computing nu(0)--> Lèvy measure has a singularity
% To avoid nu(0) in a symmetric domain, we consider N_q numero pari
y=linspace(lb,ub,N_q)'; dy=y(2)-y(1); nu_y=nu(y);
w=ones(N_q,1); w(1)=0.5; w(end)=0.5; % trapezoidal quadrature 
dx=x(2)-x(1);
dV=(V(3:end)-V(1:end-2))/(2*dx); % compute the derivative dV/dx
for i=2:length(I)-1
    I1=V_f(x,V,x(i)+y,S0,K,disc);
    I2=V(i);
    I3=(exp(y)-1)*dV(i-1);
    I(i)=sum( w.*(I1-I2-I3).*nu_y )*dy; % Numerical computation of the integral
end
end

function v=V_f(x,V,y,S0,K,disc)
% Specifies the value in the points of interest
v=zeros(size(y));
%% y<=xmin              %%%%%% BOUNDARY CONDITION (extend for x<xMIN)
%index=(y<=x(1));
%v(index)=0; %%%%%%%%%%%%%%% since it is a EU Call
%% y>= xmax             %%%%%% BOUNDARY CONDITION (extend for x>xMAX)
index=(y>=x(end));
v(index)=S0*exp(y(index))-K*disc; %%%%%% since it is a EU Call
%% xmin<y<xmax
index=find((y>x(1)).*(y<x(end)));
v(index)=interp1(x,V,y(index)); % linear interpolation
end
