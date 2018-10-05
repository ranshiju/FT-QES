function [Xi,Corr,Err]=Correlation_MPS_fitting(A,L)
Corr=zeros(L,1);

[d,D,~]=size(A);
O=zeros(d,d); O(1)=1; O(end)=0;
M0=reshape(A,[d,D*D]);
M0=reshape( permute( reshape( M0'*M0,[D,D,D,D] ),[1,3,2,4] ),[D*D,D*D] );
[vL,E]=eigs(M0.',1,'LM');
[vR,~]=eigs(M0  ,1,'LM');

Mp=reshape(A,[d,D*D]);
Mp=reshape( permute( reshape( Mp'*O*Mp,[D,D,D,D] ),[1,3,2,4] ),[D*D,D*D] );
Mp=Mp/E;

vL1=vL; vR1=vR; Z=vL'*vR; 
% Corr0=vL'*Mp*vR/Z;
for n=1:L
    vL1=Mp'*vL1; vR1=Mp*vR1;
%     Corr(n)=vL1'*vR1/Z-Corr0^2;
    Corr(n)=vL1'*vR1/Z;
end

x=(1:L)';
p=polyfit(x,log(Corr),1);Xi=abs(1/p(1));
y=exp(p(1)*x+p(2));
Err=norm(y-Corr)/norm(Corr);
keyboard
end