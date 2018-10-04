function [T,H]=TPDO_Gate_Heisenberg_Honeycomb(Jxy,Jz,hx,hz,tau)
d=2;
H=Heisenberg_TwoS(Jxy,Jz,hz/3,hx/3);

M=reshape( permute( reshape( expm(-tau*H),[d,d,d,d] ),[1,3,2,4] ),[d*d,d*d] );
[GL,lm,GR]=svd(M);

GL=reshape( permute( reshape(GL*sqrt(lm),[d,d,d*d]),[1,3,2] ),[d,d*d*d] );
GR=reshape( permute( reshape(conj(GR)*sqrt(lm),[d,d,d*d]),[1,3,2] ),[d,d*d*d] );

GL=reshape( reshape( reshape(GL,[d*d*d,d])*GL,[d*d^4,d] )*GL,[d,d*d,d*d,d*d,d] );
GR=reshape( reshape( reshape(GR,[d*d*d,d])*GR,[d*d^4,d] )*GR,[d,d*d,d*d,d*d,d] );

T=reshape( permute( reshape( reshape( permute(GL,[1,2,3,5,4]),[d*d^4*d,d*d] ) * ...
    reshape( permute(GR,[4,1,2,3,5]),[d*d,d*d^4*d] ),[d,d*d,d*d,d,d,d*d,d*d,d] ),[1,5,2,3,6,7,4,8] ),[d*d,d*d,d*d,d*d,d*d,d*d] );
end