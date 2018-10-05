function H=FullHamiltonianFiniteSize(Jxy,Jz,hx,hz,Couplings,IsPatchBh)
% IsPatchBh: 是否在边界两个site上补一半的磁场

d=2;
Nb=size(Couplings,2);
L=max(Couplings(:));

Hlocal=reshape( sptensor(Heisenberg_TwoS(Jxy,Jz,hz/2,hx/2)), [d,d,d,d]);
H=sptensor([],[],[ones(1,L)*d,ones(1,L)*d]);
Did=d^(L-2);
ID=reshape( sptensor([(1:Did)',(1:Did)'],ones(Did,1),[Did,Did]),ones(1,2*L-4)*d );
for n=1:Nb
    x1=min(Couplings(n,:)); x2=max(Couplings(n,:));
    
    % Pertmute Hbath to the right position
    V=[1,2,5:L+2,3,4,L+3:2*L];
    V=ReOrderV(V,[3:x1+1,1,x1+2:x2,2,x2+1:L,[3:x1+1,1,x1+2:x2,2,x2+1:L]+L]);
    H=H+permute(ttt(Hlocal,ID),V);
end

if(IsPatchBh)
    sx=[0,0.5;0.5,0]; sz=[0.5,0;0,-0.5];
    H1=sptensor(-hx*sx/2-hz*sz/2);
    Did=d^(L-1);
    ID=reshape( sptensor([(1:Did)',(1:Did)'],ones(Did,1),[Did,Did]),ones(1,2*L-2)*d );
    
    V=[1,3:L+1,2,L+2:2*L];
    H=H+permute(ttt(H1,ID),V);
    V=[1:L-1,2*L-1,L:2*L-2,2*L];
    H=H+permute(ttt(ID,H1),V);
end
H=reshape(H,[d^L,d^L]);
end