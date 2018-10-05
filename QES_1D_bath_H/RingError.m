function [Err,Prod]=RingError(A,v,G,d,dc,L)
tmp=reshape( reshape(G.R,[d,d^3]).'* reshape( reshape( permute( reshape(full(G.bulk),[d,d^(L-2),d,d^L]),[1,2,4,3] ),...
    [d*d^(L-2)*d^L,d] ) * reshape(G.R,[d,d^3]),[d,d^(L-2)*d^L*d^3] ),[d^2,d,d^(L-2),d^L,d^2,d] );
tmp=reshape( reshape( permute(tmp,[4,1,5,2,3,6]),[d^L*d^2*d^2,d^L] )*G.bulk,[d^L*d^2*d^2*d^L,1] );
tmp=tmp/norm(tmp)*sign(tmp(1));

A=reshape(A,[d^L,dc,dc]);
M=Get_ring_tensor(conj(A),A,permute( reshape(v.L,[dc,d*d,dc]),[2,1,3] ),permute( reshape(v.R,[dc,d*d,dc]),[2,1,3] ),d^L,d^2,dc);
M=reshape( permute(M,[1,3,2,4]),[d^L*d^2*d^2*d^L,1] );
M=M/norm(M)*sign(M(1));

Err=norm(tmp-M); Prod=tmp'*M;
end