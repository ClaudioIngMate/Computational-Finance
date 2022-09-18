function [Error,Price]=priceExtNIG_IV(x,spot,strike,rf,maturity,pmkt)

% Calibrates the model on the B&S implied volatilities

% Initialize
Price=zeros(length(strike),1);

% Characteristic exponent
char_exp=@(u) CharExpExtNIG(u,x);

% Delete duplicate
mat=unique(maturity);
for i=1:length(mat)
    index=find(maturity==mat(i));
    Price(index)=CM_price(strike(index),spot,mat(i),rf,char_exp);
end

ImpVolatilityMODEL = blsimpv(spot,strike, rf, maturity, Price);
ImpVolatilityMARKET = blsimpv(spot, strike, rf, maturity, pmkt);

Error=abs(ImpVolatilityMODEL-ImpVolatilityMARKET);

end