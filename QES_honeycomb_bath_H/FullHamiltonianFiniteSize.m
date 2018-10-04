function H=FullHamiltonianFiniteSize(Jxy,Jz,hx,hz,Couplings,PatchBh)
% PatchBh: 边界site的位置，不patch为[]

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

% Patch h/3
N=numel(PatchBh.h1);
if(N>0.1)
    sx=[0,0.5;0.5,0]; sz=[0.5,0;0,-0.5];
    H1=sptensor(-hx*sx/3-hz*sz/3);
    Did=d^(L-1);
    ID=reshape( sptensor([(1:Did)',(1:Did)'],ones(Did,1),[Did,Did]),ones(1,2*L-2)*d );
    
    for n=1:N
        V=[1,3:L+1,2,L+2:2*L];
        V=ReOrderV(V,[2:PatchBh.h1(n),1,PatchBh.h1(n)+1:L,[2:PatchBh.h1(n),1,PatchBh.h1(n)+1:L]+L]);
        H=H+permute(ttt(H1,ID),V);
    end
end

% Patch 2*h/3
N=numel(PatchBh.h2);
if(N>0.1)
    sx=[0,0.5;0.5,0]; sz=[0.5,0;0,-0.5];
    H1=sptensor(-2*hx*sx/3-2*hz*sz/3);
    Did=d^(L-1);
    ID=reshape( sptensor([(1:Did)',(1:Did)'],ones(Did,1),[Did,Did]),ones(1,2*L-2)*d );
    
    for n=1:N
        V=[1,3:L+1,2,L+2:2*L];
        V=ReOrderV(V,[2:PatchBh.h2(n),1,PatchBh.h2(n)+1:L,[2:PatchBh.h2(n),1,PatchBh.h2(n)+1:L]+L]);
        H=H+permute(ttt(H1,ID),V);
    end
end

H=reshape(H,[d^L,d^L]);
end