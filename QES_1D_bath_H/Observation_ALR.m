function [E,EL,ER,Mx,Mz]=Observation_ALR(A,AL,AR,H,d,dc,L)
% 此文件可加速
sx=[0,0.5;0.5,0]; sz=[0.5,0;0,-0.5];
E=zeros(1,L-1); Mx=zeros(1,L); Mz=Mx;
vL=reshape( eye(dc),[dc*dc,1] ); vR=vL;

UL=reshape(AL,[d^(L-1),d*dc*dc]);
UL=reshape( vL.'*reshape( permute( reshape( UL'*UL,[d,dc,dc,d,dc,dc] ),[2,5,1,4,3,6] ),[dc*dc,d*d*dc*dc] ),[d*d,dc*dc] );
UR=reshape( permute( reshape(A,[d,d^(L-1),dc*dc]),[2,1,3] ),[d^(L-1),d*dc*dc] );
UR=reshape( reshape( permute( reshape(UR'*UR,[d,dc,dc,d,dc,dc]),[2,5,1,4,3,6] ),[dc*dc*d*d,dc*dc] )*vR,[dc*dc,d*d] );
rho=reshape( permute( reshape(UL*UR,[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );
EL=trace(rho*H)/trace(rho);

UL=reshape(A,[d^(L-1),d*dc*dc]);
UL=reshape( vL.'*reshape( permute( reshape( UL'*UL,[d,dc,dc,d,dc,dc] ),[2,5,1,4,3,6] ),[dc*dc,d*d*dc*dc] ),[d*d,dc*dc] );
UR=reshape( permute( reshape(AR,[d,d^(L-1),dc*dc]),[2,1,3] ),[d^(L-1),d*dc*dc] );
UR=reshape( reshape( permute( reshape(UR'*UR,[d,dc,dc,d,dc,dc]),[2,5,1,4,3,6] ),[dc*dc*d*d,dc*dc] )*vR,[dc*dc,d*d] );
rho=reshape( permute( reshape(UL*UR,[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );
ER=trace(rho*H)/trace(rho);

A1=reshape(A,[d^L,dc*dc]);
M=reshape( permute( reshape(A1'*A1,[dc,dc,dc,dc]),[1,3,2,4] ),[dc*dc,dc*dc] );
Z=vL.'*M*vR;

for n=1:L-1
    A1=reshape(A,[d*d,d^(L-2)*dc*dc]);
    rho=A1*A1';
    E(n)=trace(rho*H);
    Mx(n)=trace(kron(eye(2),sx)*rho);
    Mz(n)=trace(kron(eye(2),sz)*rho);
    if(n==L-1)
        Mx(n+1)=trace(kron(sx,eye(2))*rho);
        Mz(n+1)=trace(kron(sz,eye(2))*rho);
    end
    
%     A1=reshape(A,[d^L,dc*dc])'*reshape( H*reshape(A,[d*d,d^(L-2)*dc*dc]),[d^L,dc*dc] );
%     E(n)=real(vL.'*reshape( permute( reshape(A1,[dc,dc,dc,dc]),[1,3,2,4] ),[dc*dc,dc*dc] )*vR/Z);
%     A=permute(reshape(A,[ones(1,L)*d,dc,dc]),[2:L,1,L+1,L+2]);
end
end