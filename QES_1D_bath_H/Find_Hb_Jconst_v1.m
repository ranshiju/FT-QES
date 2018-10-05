function [Hm,J,Error]=Find_Hb_Jconst_v1(HL,HR,d1,d2,D,tau,varargin)
Ni=numel(varargin);
if(Ni==0)
    v0=rand(1,4);
elseif(Ni==1)
    v0=varargin{1};
end
opt.TolFun=1e-12; opt.TolX=1e-12; opt.MaxFunEvals=2000; opt.maxit=5000;

[Hm]=ObtainHbath(HL,HR,d1,d2,D);
[U,LM]=eig(Hm);
LM=diag(LM)/LM(1);
LM=-log(LM)/tau;
Hm=U*diag(LM)*U';

D=min(d1*d1,d2*d2);
tmp=reshape( permute( reshape(Hm,[d1,d2,d1,d2]),[1,3,2,4] ),[d1*d1,d2*d2] );
[UL,lm,UR]=svd(tmp);
UL=reshape( UL(:,1:D)*sqrt(lm(1:D,1:D)),[d1,d1,D] ); 
UR=reshape( conj(UR(:,1:D))*sqrt(lm(1:D,1:D)),[d1,d1,D] ); 
% [U,LM]=eig(tmp);
% UL=reshape( U(:,1:D)*sqrt(LM(1:D,1:D)),[d1,d1,D] ); 
% UR=reshape( conj(U(:,1:D))*sqrt(LM(1:D,1:D)),[d1,d1,D] ); 

J.L=zeros(D,4); Error.L=zeros(D,1);
for n=1:D
    [J.L(n,:),Error.L(n)]=fminsearch(@(J0) CostF2(UL(:,:,n),J0),v0,opt);
end

J.R=zeros(D,4); Error.R=zeros(D,1);
for n=1:D
    [J.R(n,:),Error.R(n)]=fminsearch(@(J0) CostF2(UR(:,:,n),J0),v0,opt);
end
J.M=J.L.'*J.R;

H1=RestoreH(J.M);
Error.total=norm(H1-Hm);
end

function [Error]=CostF2(H,J)
% sx=[0,0.5;0.5,0]; sy=[0,1i*0.5;-1i*0.5,0]; sz=[0.5,0;0,-0.5]; I=eye(2);
% Hb1=Jb(1)*kron(sx,sx)+Jb(2)*kron(sy,sx)+Jb(2)*kron(sx,sy)+Jb(3)*kron(sy,sy)+Jb(4)*kron(sz,sz)+Jb(5)*kron(sx,sz)+Jb(5)*kron(sz,sx)+...
%     Jb(6)*kron(I,sx)+Jb(6)*kron(sx,I)+Jb(7)*kron(I,sz)+Jb(7)*kron(sz,I);

Spin={[0,0.5;0.5,0],[0,1i*0.5;-1i*0.5,0],[0.5,0;0,-0.5],eye(2)};
H1=zeros(2,2);
for n=1:4
    H1=H1+J(n)*Spin{n};
end

Error=H-H1;
Error=norm(Error(:));
end