function Price = CM_price(Strike,S0,T,r,CharExp)

% Price of a Plain Vanilla Call exploiting the Carr&Madan algorithm

% INPUT:        Strike      vector of strikes
%               S0          underlying spot price
%               T           maturity
%               r           risk-free rate
%               CharExp     characteristic exponent of the input model

% OUTPUT:       Price       vector of prices

format long

%% Discretization Parameters

Npow=15; N=2^Npow;  % discretization points
A=1000;             % domain upper boundary

%% FFT Grid (for dv and dk)

% dv
eta=A/N; v=[0:eta:A*(N-1)/N];    % integral domain grid (only > 0)
v(1)=1e-22;                      % adjust starting point near zero

% dk
lambda=2*pi/(N*eta); k=-lambda*N/2+lambda*(0:N-1);

%% FT{Z(k)} using Carr&Madan formula and inversion to find Z(k)

CharFunc=@(v) exp(T*CharExp(v)); % char function of the input model
Z_k=exp(1i*r*v*T).*(CharFunc(v-1i)-1)./(1i*v.*(1i*v+1));

w=ones(1,N); w(1)=0.5; w(end)=0.5;  % trapezoidal rule
x=w.*eta.*Z_k.*exp(1i*pi*(0:N-1));
z_k=real(fft(x)/pi);

%% EU Call price exploiting Z(k)

C=S0*(z_k+max(1-exp(k-r*T),0));     % call prices
K=S0*exp(k);                        % strikes

%%  Output: filter prices grid + interpolation

index=find( K>0.1*S0 & K<3*S0 );    % delete nonsense strikes
C=C(index); K=K(index);

Price=interp1(K,C,Strike,'spline'); % linear interpolation of price in strike

% plot(K,C)
% title( 'Option Price' );
% xlabel('Strike');
% title (' EU Call option price as function of the strike');

end
