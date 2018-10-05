function [E]=Observation_Amid(A,H,d,dc,L)
E=zeros(1,L-1);
A1=reshape(A,[d^L*dc*dc,1]);
Z=A1'*A1;

for n=1:L-1
    E(n)=real(reshape(A,[d^L*dc*dc,1])'*reshape( H*reshape(A,[d*d,d^(L-2)*dc*dc]),[d^L*dc*dc,1] )/Z);
    A=permute(reshape(A,[ones(1,L)*d,dc,dc]),[2:L,1,L+1,L+2]);
end
end