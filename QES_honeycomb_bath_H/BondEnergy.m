function Ebond=BondEnergy(v,Coupling,H,Dim)
Nb=size(Coupling,1);
N=numel(Dim);
Ebond=zeros(Nb,1); d=2;
v=reshape(v,Dim);

D1=round(prod(Dim))/(d^2);
for n=1:Nb
    x1=min(Coupling(n,1),Coupling(n,2)); x2=max(Coupling(n,1),Coupling(n,2));
    Permutation=[x1,x2,1:x1-1,x1+1:x2-1,x2+1:N];
    v1=reshape( permute(v,Permutation),[d^2,D1] );
    rho=v1*v1'; rho=rho/trace(rho);
    Ebond(n)=trace(rho*H);
end
end