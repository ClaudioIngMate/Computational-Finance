function V=CharExpMERTON(u,x)
% risk-neutral characteristic exponent


sigma=x(1);         % volatility of the diffusion component
lambda=x(2);        % intensity of Poisson process
mu=x(3);            % drift of the jump size
delta=x(4);         % standard deviation of the jump size
   
V=@(u) -sigma^2*u.^2/2+ lambda*(exp(-delta^2*u.^2/2+1i*mu*u)-1); % without drift
drift_rn=-V(-1i); % Drift Risk_neutral
V=drift_rn*1i*u+V(u);

end