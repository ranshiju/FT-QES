function [A,v,lambda,Lambda,dc2]=Redefine_v(A,v,lambda,dc,dc1,D,L)
Dh=round(sqrt(D)); d=2;
[U1,Lambda.L,V1]=svds(reshape(v.L,[dc1*Dh,Dh*dc1]),dc);
[U2,Lambda.R,V2]=svds(reshape(v.R,[dc1*Dh,Dh*dc1]),dc);

dc2=min(numel(diag(Lambda.L)),numel(diag(Lambda.R)));
U1=U1(:,1:dc2); V1=V1(:,1:dc2); Lambda.L=Lambda.L(1:dc2,1:dc2);
U2=U2(:,1:dc2); V2=V2(:,1:dc2); Lambda.R=Lambda.R(1:dc2,1:dc2);

v.L=reshape( reshape((V1*sqrt(Lambda.L))',[dc2*Dh,dc1]) * reshape(U1*sqrt(Lambda.L),[dc1,Dh*dc2]),[dc2*D*dc2,1] );
v.R=reshape( reshape((V2*sqrt(Lambda.R))',[dc2*Dh,dc1]) * reshape(U2*sqrt(Lambda.R),[dc1,Dh*dc2]),[dc2*D*dc2,1] );

if(dc2<dc1)
    A=reshape(A,[d^L,dc1,dc1]); A=reshape(A(:,1:dc2,1:dc2),[d^L*dc2*dc2,1]);
    lambda=lambda(1:dc2);
elseif(dc2>dc1)
    [A]=Initial_tensor(L,dc2,d,2); lambda=ones(dc2,1);
end

Lambda.L=diag(Lambda.L); Lambda.R=diag(Lambda.R);
end