function [Xi]=Correlation_iDMRG(AL,AR,L,dc)
d=2; Xi=zeros(1,2);
M=reshape(AL,[d^L,dc*dc]);
[D]=eigs(reshape( permute( reshape( M'*M,[dc,dc,dc,dc] ),[2,4,1,3] ),[dc*dc,dc*dc] ),2,'LM');
D=abs(D(1)/D(2)); Xi(1)=1/log(D);

M=reshape(AR,[d^L,dc*dc]);
[D]=eigs(reshape( permute( reshape( M'*M,[dc,dc,dc,dc] ),[1,3,2,4] ),[dc*dc,dc*dc] ),2,'LM');
D=abs(D(1)/D(2)); Xi(2)=1/log(D);

Xi=Xi*L;
end