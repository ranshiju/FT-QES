function [Mlocal,GSlm,GSEnt]=ObservationGS(v,Nspin,d,D)
S=size(v); 
% Mstate=S(2);
sx=[0,0.5;0.5,0]; sy=[0,0.5i;-0.5i,0]; sz=[0.5,0;0,-0.5];
% Ob.Mtot=zeros(3,Mstate);
Mlocal=zeros(3,Nspin-2);

% Now only calculate the ground state
% for m=1:Mstate
Norm=v(:,1)'*v(:,1);
v1=reshape(v(:,1),[D,ones(1,Nspin-2)*d,D])/Norm;
D1=D^2*d^(Nspin-3);
for n=2:Nspin-1
%         Ob.Mtot(1,m)=Ob.Mtot(1,m)+conj(vm)*reshape(sx*v1,[d^(Nspin),1])/Norm;
%         Ob.Mtot(2,m)=Ob.Mtot(2,m)+conj(vm)*reshape(sy*v1,[d^(Nspin),1])/Norm;
%     if(m==1)
        vv=reshape( permute(v1,[n,1:n-1,n+1:Nspin]),[d,D1] );
        vv=vv*vv';
        Mlocal(1,n-1)=trace(vv*sx);
        Mlocal(2,n-1)=trace(vv*sy);
        Mlocal(3,n-1)=trace(vv*sz);
%             Ob.Mtot(3,m)=Ob.Mtot(3,m)+Ob.Mzlocal(n);
%     else
%             Ob.Mtot(3,m)=Ob.Mtot(3,m)+conj(vm)*reshape(sz*v1,[d^(Nspin),1])/Norm;
%     end
end
% end

if(rem(Nspin,2)==0 && Nspin-2>0.1)
    Dh=D*d^round((Nspin-2)/2);
    GSlm=svd(reshape(v(:,1),[Dh,Dh]));
    GSlm=GSlm/norm(GSlm);
    D=GSlm.^2;
    GSEnt=-D'*log(D);
end
end