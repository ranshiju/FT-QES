function [E,Econnect,vL,vR,lambda]=Observe_Energy(A,H,L,dc,varargin)
% Input: (d,d,...,d,dc,dc)
Ni=numel(varargin);

d=2; E=zeros(1,L-1);
A1=reshape(A,[d^L,dc*dc]);
M=reshape( permute( reshape(A1'*A1,[dc,dc,dc,dc]),[1,3,2,4] ),[dc*dc,dc*dc] );

if(Ni==0)
    opt.v0=reshape(eye(dc),dc*dc,1); opt.isreal=0;
    [vL,~]=eigs(M.',1,'LM',opt); [vR,~]=eigs(M,1,'LM',opt);
elseif(Ni==2)
    vL=reshape(varargin{1},[dc*dc,1]); 
    vR=reshape(varargin{2},[dc*dc,1]);
else
    error('Bad input in function Observe_Energy.m! \n')
end
[~,~,lambda]=TruncationFromRM(reshape(vL,[dc,dc]),reshape(vR,[dc,dc]),ones(dc,1),0);

% vL=reshape(eye(dc),dc*dc,1); vR=vL;

Z=vL.'*M*vR;
d=2;

AL=reshape(A,[d^(L-1),d*dc*dc]);
AL=reshape( vL.'*reshape( permute( reshape( AL'*AL,[d,dc,dc,d,dc,dc] ),[2,5,1,4,3,6] ),[dc*dc,d*d*dc*dc] ),[d*d,dc*dc] );
AR=reshape( permute( reshape(A,[d,d^(L-1),dc*dc]),[2,1,3] ),[d^(L-1),d*dc*dc] );
AR=reshape( reshape( permute( reshape(AR'*AR,[d,dc,dc,d,dc,dc]),[2,5,1,4,3,6] ),[dc*dc*d*d,dc*dc] )*vR,[dc*dc,d*d] );
rho=reshape( permute( reshape(AL*AR,[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );
Econnect=trace(rho*H)/trace(rho);

for n=1:L-1
    A1=reshape(A,[d^L,dc*dc])'*reshape( H*reshape(A,[d*d,d^(L-2)*dc*dc]),[d^L,dc*dc] );
    E(n)=real(vL.'*reshape( permute( reshape(A1,[dc,dc,dc,dc]),[1,3,2,4] ),[dc*dc,dc*dc] )*vR/Z);
    A=permute(reshape(A,[ones(1,L)*d,dc,dc]),[2:L,1,L+1,L+2]);
end


end