function [Mlocal,Eb,Z]=ObservationFT(rho,Hbulk,Hbath,Nspin,d,D)
sx=[0,0.5;0.5,0]; sy=[0,0.5i;-0.5i,0]; sz=[0.5,0;0,-0.5];
Mlocal=zeros(3,Nspin-2);
Eb=zeros(1,Nspin-1);

D1=D^2*d^(Nspin-3);
Z=trace(rho);
rho=reshape(rho,[D,ones(1,Nspin-2)*d,D,D,ones(1,Nspin-2)*d,D]);
Id=reshape( eye(D1,D1),[D1*D1,1] );

for n=2:Nspin-1
        rho1=reshape( reshape( permute(rho,[n,n+Nspin,1:n-1,n+1:Nspin,Nspin+1:Nspin+n-1,Nspin+n+1:Nspin*2]),[d*d,D1*D1] )*Id,[d,d] );
        Mlocal(1,n-1)=trace(rho1*sx)/Z;
        Mlocal(2,n-1)=trace(rho1*sy)/Z;
        Mlocal(3,n-1)=trace(rho1*sz)/Z;
end

for n=1:Nspin-1
    if(n==1 || n==Nspin-1)
        D1=D*d^(Nspin-3);
        Id=reshape( eye(D1,D1),[D1*D1,1] );
        rho1=reshape( reshape( permute(rho,[n,n+1,n+Nspin,n+1+Nspin,1:n-1,n+2:Nspin,Nspin+1:Nspin+n-1,Nspin+n+2:Nspin*2]),[d*D*d*D,D1*D1] )*Id,[d*D,d*D] );
        Eb(1,n)=trace(rho1*Hbath{(n==1)+(n==Nspin-1)*2})/Z;
    else
        D1=D^2*d^(Nspin-4);
        Id=reshape( eye(D1,D1),[D1*D1,1] );
        rho1=reshape( reshape( permute(rho,[n,n+1,n+Nspin,n+1+Nspin,1:n-1,n+2:Nspin,Nspin+1:Nspin+n-1,Nspin+n+2:Nspin*2]),[d*d*d*d,D1*D1] )*Id,[d*d,d*d] );
        Eb(1,n)=trace(rho1*Hbulk)/Z;
    end
end
end