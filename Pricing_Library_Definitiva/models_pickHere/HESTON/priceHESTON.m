function [Error,Price]=priceHESTON(x,S0,strike,r,maturity,pmkt)

% Initialize
Price=zeros(length(strike),1);

% Delete duplicate
mat=unique(maturity);
for i=1:length(mat)
    index=find(maturity==mat(i));
    Price(index)=CM_priceHESTON(x,S0,strike(index),r,mat(i));
end

Error=abs(pmkt-Price)./pmkt;

end