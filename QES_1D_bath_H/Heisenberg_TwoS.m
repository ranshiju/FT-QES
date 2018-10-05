function H=Heisenberg_TwoS(jxy,jz,hz,hx)
su=[0,1;0,0];
sd=[0,0;1,0];
sz=[0.5,0;0,-0.5];
sx=[0,0.5;0.5,0];
% sy=[0,0.5*1i;-0.5*1i,0];

I=eye(2,2);
% sx=[0,1;1,0];
% d=2; H=reshape( jxy*(kron(su,sd)+kron(sd,su))/2+jz*kron(sz,sz)-h/2*(kron(sz,I)+kron(I,sz)),[d,d,d,d] );
H=jxy*(kron(su,sd)+kron(sd,su))/2+jz*kron(sz,sz)-hz*(kron(sz,I)+kron(I,sz))-hx*(kron(sx,I)+kron(I,sx));
% H=jxy*(kron(sx,sx)+kron(sy,y))+jz*kron(sz,sz)-hz*(kron(sz,I)+kron(I,sz))-hx*(kron(sx,I)+kron(I,sx));
end