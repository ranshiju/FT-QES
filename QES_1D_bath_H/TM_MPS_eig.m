function [UL,UR]=TM_MPS_eig(A,d,dc,L)
A1=reshape(A,[d^L,dc*dc]);
M=reshape( permute( reshape( A1'*A1,[dc,dc,dc,dc] ),[1,3,2,4] ),[dc*dc,dc*dc] );
[UL,~]=eigs(M.',1,'LM'); UL=reshape(UL,[dc,dc]);
[UR,~]=eigs(M  ,1,'LM'); UR=reshape(UR,[dc,dc]);
end