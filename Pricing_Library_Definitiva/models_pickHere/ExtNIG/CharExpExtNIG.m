function V=CharExpExtNIG(u,x)
% risk-neutral characteristic exponent

% Model parameters
sigma=x(1);         % volatility of the diffusion component
theta=x(2);         % drift of the diffusion component
k=x(3);             % variance of the subordinator
sigmaGBM=x(4);      % volatility of the Pure B.M. (No time change)

V=@(u) -(sigmaGBM^2/2)*u.^2 +1/k*(1-sqrt(1+u.^2*sigma.^2*k-2*1i*theta*k*u)); % without drift
drift_rn=-V(-1i); % Drift Risk_neutral
V=drift_rn*1i*u+V(u);

end