function Hm=Obtain_Hfb(HL,HR,d1,d2,tau)
[Hm]=ObtainHbath(HL,HR,d1,d2,4);
[U,LM]=eig(Hm);
LM=diag(LM)/LM(1);
D=numel(LM);
LM=(ones(size(LM))-LM)/tau;
Hm=U*diag(LM)*U';
Hm=real(Hm-trace(Hm)/D*eye(size(Hm)));
end