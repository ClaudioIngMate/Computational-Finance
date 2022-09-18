function [Error,Price]=priceMERTON(x,S0,strike,r,maturity,pmkt)

% If I want only to compute the C&M price, without calibrating anything
if nargin==5
    pmkt=nan;
end

% Initialize
Price=zeros(length(strike),1);

% Characteristic Exponent
char_exp=@(u) CharExpMERTON(u,x);

% Delete duplicate
mat=unique(maturity);
for i=1:length(mat)
    index=find(maturity==mat(i));
    Price(index)=CM_price(strike(index),S0,mat(i),r,char_exp);
end

Error=abs(pmkt-Price)./pmkt;

end