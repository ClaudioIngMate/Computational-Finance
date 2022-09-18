function V=CharExpKOU(u,x)
% risk-neutral characteristic exponent

% Model parameters
sigma=x(1);             % volatility of the diffusion component
lambda=x(2);            % intensity of Poisson Process
p=x(3);                 % probability of positive jumps
lambdaplus=x(4);        % law of positive jumps
lambdaminus=x(5);       % law of negative jumps    

V=@(u) -sigma^2*u.^2/2+...
    1i*lambda*u.*(p./(lambdaplus-1i*u)-(1-p)./(lambdaminus+1i*u)); % without drift
drift_rn=-V(-1i); % Drift Risk_neutral
V=drift_rn*1i*u+V(u);

end