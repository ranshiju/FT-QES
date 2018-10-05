function [A,S]=Stablizer_entanglement(A,v,G,L,dc,Way)
d=2;
if(Way==1)
    D=d^2;
elseif(Way==2)
    D=d^4;
end

M=TM_Update_v(A,G.L,G.R,G.bulk,L,d,dc,Way);
v.L=M.'*v.L; v.L=v.L/norm(v.L);
v.R=M*v.R;   v.R=v.R/norm(v.R);

v.L=reshape( permute( reshape(v.L,[dc,D,dc]),[1,3,2] ),[dc*dc,D] );
v.R=reshape( permute( reshape(v.R,[dc,D,dc]),[2,1,3] ),[D,dc*dc] );
M=reshape( permute( reshape( v.L*v.R,[dc,dc,dc,dc] ),[1,3,2,4] ),[dc*dc,dc*dc] );

% opt.isreal=0; 
[M,~]=eigs(M,1);
% M=reshape(eye(dc),[1,dc*dc])*M;

[U,S,V]=svd(reshape(M,[dc,dc]));
U=U*sqrt(S); V=V*sqrt(S);

A=reshape( V'*reshape( permute( reshape(A,[d^L,dc,dc]),[2,1,3] ),[dc,d^L*dc] ),[dc*d^L,dc] )*U;
A=reshape( permute( reshape(A,[dc,d^L,dc]),[2,1,3] ),[d^L*dc*dc,1] );
S=diag(S);
end