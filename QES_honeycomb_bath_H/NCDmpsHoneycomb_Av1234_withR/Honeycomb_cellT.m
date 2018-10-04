function [T,H2,U]=Honeycomb_cellT(Jxy,Jz,hx,hz,tau,Order)
d=2;

H2=Heisenberg_TwoS(Jxy,Jz,hz/3,hx/3);
% G=expm( -tau*H2 );
G=eye(d*d)-tau*H2;
T=reshape( permute( reshape(G,[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );
[U{1},LM,U{2}]=svd(T);
U{1}=U{1}*sqrt(LM);
U{2}=conj(U{2})*sqrt(LM);

UL=permute( reshape(U{1},[d,d,d*d]),[1,3,2] );
UL=reshape( reshape(UL,[d*d*d,d]) * reshape(UL,[d,d*d*d]),[d,d*d,d*d,d] );
UR=permute( reshape(U{2},[d,d,d*d]),[1,3,2] );
UR=reshape( reshape(UR,[d*d*d,d]) * reshape(UR,[d,d*d*d]),[d,d*d,d*d,d] );
G=reshape(G,[d,d,d,d]);

if(Order==2)
    T=ncon({UL,UR,G,G},{[1,-3,-4,2],[3,-5,-6,4],[-1,-2,1,3],[2,4,-7,-8]});
elseif(Order==1)
    T=ncon({UL,UR,G},{[1,-3,-4,-7],[2,-5,-6,-8],[-1,-2,1,2]});
end
T=reshape(T,ones(1,6)*d*d);
end