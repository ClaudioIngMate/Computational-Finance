function [x,h,w,H] = kernel(ngrid,xmin,xmax,params,dt,flag)

% INPUT         ngrid           number of grid points
%               [xmin, xmax]    interval where we compute the kernel 
%               params          Model parameters
%               flag            0 for backward (default) , 1 for forward problem

if nargin==4
    flag=0;         % backward problem by default
end

%% Create the grid

N = ngrid/2;
dx = (xmax-xmin)/ngrid;
x = dx*(-N:N-1);
dw = 2*pi/(xmax-xmin);
w = dw*(-N:N-1);

%% NIG characteristic exponent
sigma=params(1); theta=params(2); k=params(3);
% Characteristic exponent without drift
V=@(u) 1/k*(1-sqrt(1+u.^2*sigma.^2*k-2*1i*theta*k*u));

% %% Extended NIG characteristic exponent
% sigma=params(1); theta=params(2); k=params(3); sigmaGBM=params(4);
% % Characteristic exponent without drift 
% V=@(u) -(sigmaGBM^2/2)*u.^2 +1/k*(1-sqrt(1+u.^2*sigma.^2*k-2*1i*theta*k*u));

% %% VG characteristic exponent
% sigma=params(1); theta=params(2); k=params(3);
% % Characteristic exponent without drift 
% V=@(u) -(1/k)*log(1+u.^2*sigma.^2*k/2-1i*theta*k*u);

% %% Extended VG characteristic exponent
% sigma=params(1); theta=params(2); k=params(3); sigmaGBM=params(4);
% % Characteristic exponent without drift
% V=@(u) -(sigmaGBM^2/2)*u.^2 -1/k*log(1+u.^2*sigma.^2*k/2-1i*theta*k*u);

% %% KOU characteristic exponent
% sigma=params(1); lambda=params(2); p=params(3);
% lambdaplus=params(4); lambdaminus=params(5);
% % Characteristic exponent without drift
% V=@(u) -sigma^2*u.^2/2+ 1i*lambda*u.*(p./(lambdaplus-1i*u)-(1-p)./(lambdaminus+1i*u));

% %% MERTON characteristic exponent
% sigma=params(1); lambda=params(2); mu=params(3); delta=params(4);
% % Characteristic exponent without drift
% V=@(u) -sigma^2*u.^2/2+ lambda*(exp(-delta^2*u.^2/2+1i*mu*u)-1);

%% Characteristic function computation (NOT use with HESTON model)

drift_rn=-V(-1i);           % risk - neutral drift
V=drift_rn*1i*w+V(w);       % risk - neutral char. exponent
H=exp(dt*V);                % characteristic function


% %% HESTON characteristic function 
% 
% theta=params(1); csi=params(2); rho=params(3); eta_=params(4); V0=params(5);
% u=w; T=dt;
% alfa=-.5*(u.*u+u*1i); beta=csi-rho*theta*u*1i; epsilon2=theta*theta; gamma=.5*epsilon2;
% D=sqrt(beta.*beta-4.0*alfa.*gamma); bD=beta-D; eDt=exp(-D*T); G=bD./(beta+D);
% B=(bD./epsilon2).*((1.0-eDt)./(1.0-G.*eDt));
% psi=(G.*eDt-1.0)./(G-1.0); A_=((csi*eta_)/(epsilon2))*(bD*T-2.0*log(psi));
% y = A_ + B*V0;
% H = exp(y);

%%
if flag==0
    H=conj(H);              % conjugate characteristic function
end

%% Kernel computation

h = real(fftshift(fft(ifftshift(H))))/(xmax-xmin); % kernel

% figure
% plot(x,h)
% figure
% plot(w,H)

end