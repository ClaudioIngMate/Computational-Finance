function [Error,Price]=priceBATES_IV(x,S0,strike,r,maturity,pmkt)

% Calibrates the model on the B&S implied volatilities

% Initialize
Price=zeros(length(strike),1);

% Delete duplicate
mat=unique(maturity);
for i=1:length(mat)
    index=find(maturity==mat(i));
    Price(index)=CM_priceBATES(x,S0,strike(index),r,mat(i));
end

ImpVolatilityMODEL = blsimpv(S0,strike, r, maturity, Price.*(Price>0));
ImpVolatilityMARKET = blsimpv(S0, strike, r, maturity, pmkt.*(pmkt>0));

Error=abs(ImpVolatilityMODEL-ImpVolatilityMARKET);

end