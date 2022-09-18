function [Error,Price]=priceBATES_IV(x,spot,strike,rf,maturity,pmkt)

% Initialize
Price=zeros(length(strike),1);



% Delete duplicate
mat=unique(maturity);
for i=1:length(mat)
    index=find(maturity==mat(i));
    Price(index)=CM_priceBATES(x,spot,strike(index),rf,mat(i));
end

ImpVolatilityMODEL = blsimpv(spot,...
    strike, rf, maturity, Price);
ImpVolatilityMARKET = blsimpv...
    (spot, strike, rf, maturity, pmkt);
Error=abs(ImpVolatilityMODEL-ImpVolatilityMARKET);

end