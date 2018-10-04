function [Hx,Hy,Hz]=HamiltonianKitaHesen(Jh,Jk,hx,hz)
sx=[0,0.5;0.5,0]; sy=[0,0.5*1i;-0.5*1i,0]; sz=[0.5,0;0,-0.5];
H=Heisenberg_TwoS(Jh,Jh,hz/3,hx/3);

Hx=H+kron(sx,sx)*Jk;
Hy=H+real(kron(sy,sy))*Jk;
Hz=H+kron(sz,sz)*Jk;
end