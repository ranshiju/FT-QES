function Ebond=BondEnergyFL(v,Coupling,H,N)
Nb=size(Coupling,1);
Ebond=zeros(Nb,1); d=2;
v=reshape(v,ones(1,N)*d);

for n=1:Nb
    x1=min(Coupling(n,1),Coupling(n,2)); x2=max(Coupling(n,1),Coupling(n,2));
    Permutation=[x1,x2,1:x1-1,x1+1:x2-1,x2+1:N];
    v1=reshape( permute(v,Permutation),[d^2,d^(N-2)] );
    rho=v1*v1'; rho=rho/trace(rho);
    Ebond(n)=trace(rho*H);
end
end