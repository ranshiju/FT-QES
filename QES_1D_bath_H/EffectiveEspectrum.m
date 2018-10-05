function [Erho,Es,Xiv]=EffectiveEspectrum(v,G,H,L,d,dc,para)
Es=zeros(1,para.Nstate);

Hhandle=@(A) Update_A_handle(A,v,G,para);
[Phi,Erho]=eigs(Hhandle,d^L*dc*dc,para.Nstate,'LM');
Erho=diag(Erho).';
tmp=abs(real(Erho(1))/real(Erho(2))); Xiv=1/log(tmp);

for n=1:para.Nstate
    tmp=Observation_Amid(Phi(:,n),H,d,dc,L);
    Es(n)=sum(tmp)/(L-1);
end
Es=sort(Es,'ascend');
end