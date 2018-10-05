function [Mlocal,Eb,Z]=ObservationFTFL(rho,Hbulk,Nspin,d)
% For the observation of a finite-size chain without entanglement bath
% 注意：观测Eb时，将磁场平分到两个site上，会导致边界两个site的磁场较小一半，故Eb的求和比总能量小2*h/2*mz_local

sx=[0,0.5;0.5,0]; sy=[0,0.5i;-0.5i,0]; sz=[0.5,0;0,-0.5];
Mlocal=zeros(3,Nspin);
Eb=zeros(1,Nspin-1);

Z=trace(rho);
rho=reshape(rho,ones(1,2*Nspin)*d);
D1=d^(Nspin-1); 
Id=reshape( eye(D1,D1),[D1*D1,1] );

for n=1:Nspin
    rho1=reshape( reshape( permute(rho,[n,n+Nspin,1:n-1,n+1:Nspin,Nspin+1:Nspin+n-1,Nspin+n+1:Nspin*2]),[d*d,D1*D1] )*Id,[d,d] );
    Mlocal(1,n)=trace(rho1*sx)/Z;
    Mlocal(2,n)=trace(rho1*sy)/Z;
    Mlocal(3,n)=trace(rho1*sz)/Z;
end

D1=d^(Nspin-2);
Id=reshape( eye(D1,D1),[D1*D1,1] );
for n=1:Nspin-1
    rho1=reshape( reshape( permute(rho,[n,n+1,n+Nspin,n+1+Nspin,1:n-1,n+2:Nspin,Nspin+1:Nspin+n-1,Nspin+n+2:Nspin*2]),[d*d*d*d,D1*D1] )*Id,[d*d,d*d] );
    Eb(1,n)=trace(rho1*Hbulk)/Z;
end
end