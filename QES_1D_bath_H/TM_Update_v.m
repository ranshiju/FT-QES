function [M,A,v,dc2]=TM_Update_v(A,v,GL,GR,Gbulk,L,d,dc,Otrotter)
% Input: A (d,d,...,d,dc,dc)

dc2=dc;
% [A,lambda,dc2]=RedefineA(A,d,dc,dc1,L,lambda,1);

A1=Evolve_H_Bulk(A,Gbulk,L,d,dc2,1);
A1=reshape(permute(reshape(A1,[ones(1,L)*d,dc2,dc2]),[2:L,1,L+1,L+2]),[d^(L-2),d*d*dc2*dc2]);
if(Otrotter==2)
    M=permute( reshape( A1'*A1,[d,d,dc2,dc2,d,d,dc2,dc2] ),[2,6,3,7,4,8,1,5] );
elseif(Otrotter==1)
    tmp=reshape( permute( reshape(A,[ones(1,L)*d,dc2,dc2]),[2:L,1,L+1,L+2] ),[d^(L-2),d*d*dc2*dc2]);
    M=permute( reshape( A1'*tmp,[d,d,dc2,dc2,d,d,dc2,dc2] ),[2,6,3,7,4,8,1,5] );
end
GL=reshape( permute( reshape(GL,[d,d*d,d]),[1,3,2] ),[d*d,d*d] );
GR=reshape( permute( reshape(GR,[d,d*d,d]),[2,1,3] ),[d*d,d*d] );
D=d*d;

M=GR*reshape(M,[d*d,dc2*dc2*dc2*dc2*d*d]);
M=reshape( reshape(M,[D*dc2*dc2*dc2*dc2,d*d])*GL,[D,dc2,dc2,dc2,dc2,D] );
M=reshape( permute(M,[2,1,3,4,6,5]),[dc2*D*dc2,dc2*D*dc2] );

if(dc2<dc)
    v.L=reshape(v.L,[dc,D,dc]);  v.L=reshape( v.L(1:dc2,:,1:dc2),[dc2*D*dc2,1] );
    v.R=reshape(v.R,[dc,D,dc]);  v.R=reshape( v.R(1:dc2,:,1:dc2),[dc2*D*dc2,1] );
elseif(dc2>dc)
    [~,v]=Initial_tensor(L,dc2,d,Way);
end
end