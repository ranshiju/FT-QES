function [Hb,Jb,Error]=Find_Hb_Jconst(HL,HR,d1,d2,D,tau)
[Hb]=ObtainHbath(HL,HR,d1,d2,D);
[U,LM]=eig(Hb);
LM=diag(LM)/LM(1);
LM=-log(LM)/tau;
Hb=U*diag(LM)*U';

if(d1==2 && d2==2)
    [Jb,Error]=fminsearch(@(Jb0) CostF2(Hb,Jb0),randn(16,1));
end

end

function [Error]=CostF2(Hb,Jb)
% sx=[0,0.5;0.5,0]; sy=[0,1i*0.5;-1i*0.5,0]; sz=[0.5,0;0,-0.5]; I=eye(2);
% Hb1=Jb(1)*kron(sx,sx)+Jb(2)*kron(sy,sx)+Jb(2)*kron(sx,sy)+Jb(3)*kron(sy,sy)+Jb(4)*kron(sz,sz)+Jb(5)*kron(sx,sz)+Jb(5)*kron(sz,sx)+...
%     Jb(6)*kron(I,sx)+Jb(6)*kron(sx,I)+Jb(7)*kron(I,sz)+Jb(7)*kron(sz,I);

Spin={[0,0.5;0.5,0],[0,1i*0.5;-1i*0.5,0],[0.5,0;0,-0.5],eye(2)};
Jb=reshape(Jb,[4,4]);
Hb1=zeros(4,4);
for n1=1:4
    for n2=1:4
        Hb1=Hb1+Jb(n1,n2)*kron(Spin{n1},Spin{n2});
    end
end
Jb=reshape(Jb,[16,1]);

Error=Hb-Hb1;
Error=norm(Error(:));
end