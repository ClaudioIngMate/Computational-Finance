function [Error,Price]=priceVG_IV(x,S0,strike,r,maturity,pmkt)

% Calibrates the model on the B&S implied volatilities

% Initialize
Price=zeros(length(strike),1);

% Characteristic exponent
char_exp=@(u) CharExpVG(u,x);

% Delete duplicate
mat=unique(maturity);
for i=1:length(mat)
    index=find(maturity==mat(i));
    Price(index)=CM_price(strike(index),S0,mat(i),r,char_exp);
end

ImpVolatilityMODEL = blsimpv(S0,strike, r, maturity, Price);
ImpVolatilityMARKET = blsimpv(S0, strike, r, maturity, pmkt);

Error=abs(ImpVolatilityMODEL-ImpVolatilityMARKET);

end