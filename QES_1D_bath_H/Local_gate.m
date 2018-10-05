function [G,H,Hbulk]=Local_gate(Jxy,Jz,hx,hz,tau,L,Issparse,Otrotter)
d=2;

H=Heisenberg_TwoS(Jxy,Jz,hz/2,hx/2);
G.G2=eye(size(H))-tau*H;

M=reshape( permute( reshape(G.G2,[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );
[G.L,D,G.R]=svd(M);
G.L=reshape( permute( reshape(G.L*sqrt(D),[d,d,d*d]),[1,3,2] ),[d*d*d,d] );
G.R=reshape( permute( reshape(G.R*sqrt(D),[d,d,d*d]),[1,3,2] ),[d*d*d,d] );

if(Issparse)
    Hbulk=Hamiltonian_Chain_sparse(H,L);
    G.bulk=sparse(eye(d^L))-tau*Hbulk/(1+(Otrotter==2));
else
    Hbulk=[];
end
end