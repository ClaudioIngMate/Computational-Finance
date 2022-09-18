function [S,v] = CONV(S0, K, r, T, Ndates, N, Barrier_L, Barrier_U, param)

% INPUT     S0 (spot price), K (strike), Ndates ( # monitoring dates)
%           N (grid length), Barrier_L (lower) , Barrier_U (upper)
%           param (vector with model parameters)

%OUTPUT
%       S: grid of spot prices at maturity (ST)
%       v: grid of values of the Down and Out Call option      

dt=T/Ndates;
b=2.5; % truncation domain. Means Smin/max=S0*exp( -/+ b)-> [0.08, 12.18] if S0=1;
[x,~,~,H] = kernel(N,-b,b,param,dt,0); % H contains the conjugate char. function
S = S0*exp(x);

v = max(S-K,0).*(S>Barrier_L).*(S<Barrier_U); % Payoff at maturity

H=ifftshift(H); % we move to Matlab view
for j = 1:Ndates
    %v(S<=Barrier_L) = 0;  % if in t=0 Barrier has NO effect
    %v(S>=Barrier_U) = 0;  % if in t=0 Barrier has NO effect
    v=real(fftshift(fft(ifft(ifftshift(v)).*H)))*exp(-r*dt);
    v(S<=Barrier_L) = 0;   % if in t=0 Barrier has effect
    v(S>=Barrier_U) = 0;   % if in t=0 Barrier has effect
end

index=find( (S>0.1*S0).*(S<3*S0));
S=S(index); v=v(index);
figure
plot(S,v,'*', 'markersize',0.8);
xlabel('S0'); title('Knock & Out Call price'); grid on;

end