function [ent]=Entanglement(lambda)
S=size(lambda);
if(S(1)>1.1 && S(2)>1.1)
    warning('Input should be a column vector instead of matrix')
elseif(S(1)<S(2)-0.1)
    lambda=lambda.'; % Make sure it is a column vector
end

lambda=lambda/norm(lambda);
Eps=1e-20;
dc=find(lambda>Eps,1,'last');
lambda=lambda(1:dc);
ent=-(lambda.^2)'*log(lambda.^2);
end