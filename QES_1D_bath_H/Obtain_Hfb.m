function Hm=Obtain_Hfb(HL,HR,d1,d2,tau)
[Hm]=ObtainHbath(HL,HR,d1,d2,4);
[U,LM]=eig(Hm);
LM=diag(LM)/LM(1);
LM=-log(LM)/tau;
Hm=U*diag(LM)*U';
Hm = Hm - eye(size(Hm)) / size(Hm, 1) * trace(Hm);
end