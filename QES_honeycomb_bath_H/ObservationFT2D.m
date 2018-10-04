function [Mlocal,Eb,Z]=ObservationFT2D(rho,Hbulk,PosP,PosB,Ih2,Dim,d,D)
Np=numel(PosP);
Nb=numel(PosB);
Nc=size(Ih2,1);

sx=[0,0.5;0.5,0]; sy=[0,0.5i;-0.5i,0]; sz=[0.5,0;0,-0.5];
Mlocal=zeros(3,Np);

D1=D^Nb*d^(Np-1);
Z=trace(rho);
rho=reshape(rho,[Dim,Dim])/Z;
Id=reshape( eye(D1,D1),[D1*D1,1] );

for n=1:Np
    Bond=[PosP(n),PosP(n)+Np+Nb,1:PosP(n)-1,PosP(n)+1:Np+Nb,[1:PosP(n)-1,PosP(n)+1:Np+Nb]+Np+Nb];
    rho1=reshape( reshape( permute(rho,Bond),[d*d,D1*D1] )*Id,[d,d] );
    Mlocal(1,n-1)=trace(rho1*sx);
    Mlocal(2,n-1)=trace(rho1*sy);
    Mlocal(3,n-1)=trace(rho1*sz);
end

Eb=zeros(1,Nc);
D1=D^Nb*d^(Np-2);
for n=1:Nc
    x1=Ih2(n,1); x2=Ih2(n,2);
    Id=reshape( eye(D1,D1),[D1*D1,1] );
    Bond=[x1,x2,[x1,x2]+Np+Nb,1:x1-1,x1+1:x2-1,x2+1:Np+Nb,[1:x1-1,x1+1:x2-1,x2+1:Np+Nb]+Np+Nb];
    rho1=reshape( reshape( permute(rho,Bond),[d*d*d*d,D1*D1] )*Id,[d*d,d*d] );
    Eb(1,n)=trace(rho1*Hbulk);
end
end