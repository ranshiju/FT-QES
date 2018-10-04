function [Mlocal,GSlm,GSEnt]=ObservationGS(v,Dim,PosP,PosB)
% Only calculate the ground state
% No calculating entanglement

sx=[0,0.5;0.5,0]; sy=[0,0.5i;-0.5i,0]; sz=[0.5,0;0,-0.5];
Lp=numel(PosP);
Lb=numel(PosB);
L=Lp+Lb;
d=Dim(PosP(1));
D=Dim(PosB(1));

Norm=v';
v=reshape(v,Dim)/Norm;

Mlocal=zeros(3,Lp);
D1=D^Lb*d^(Lp-1);
for n=1:Lp
    vv=reshape( permute(v,[PosP(n),1:PosP(n)-1,PosP(n)+1:L]),[d,D1] );
    vv=vv*vv';
    Mlocal(1,n-1)=trace(vv*sx);
    Mlocal(2,n-1)=trace(vv*sy);
    Mlocal(3,n-1)=trace(vv*sz);
end

GSlm=[]; GSEnt=[];
% if(rem(L,2)==0 && L-2>0.1)
%     Dh=D*d^round((Nspin-2)/2);
%     GSlm=svd(reshape(v(:,1),[Dh,Dh]));
%     GSlm=GSlm/norm(GSlm);
%     D=GSlm.^2;
%     GSEnt=-D'*log(D);
% end
end