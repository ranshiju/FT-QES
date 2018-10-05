function [Xi]=Correlation_MPS(A,d,dc)
M=reshape(A,[d,dc*dc]);
[D]=eigs(reshape( permute( reshape( M'*M,[dc,dc,dc,dc] ),[1,3,2,4] ),[dc*dc,dc*dc] ),2,'LM');
D=abs(D(1)/D(2)); Xi=1/log(D);
end