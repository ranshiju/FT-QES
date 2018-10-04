function [v4,A,R,LM,Norm]=Update_v4(T,A,v1,v2,v3)
[dp,D,~,~,~,~]=size(T);
[~,d1,d2,d3,d4]=size(A);

% [A,R]=qr(reshape(A,[dp*d1*d2*d3,d4]),0);
% A=reshape(A,[dp,d1,d2,d3,d4]); R=R/norm(R(:));
% v4=reshape( reshape(T,[dp*D^4,dp]) * reshape(A,[dp,d1*d2*d3*d4]),[dp,D,D,D,D,d1,d2,d3,d4] );

[A,LM,V]=svd(reshape(A,[dp*d1*d2*d3,d4]),0); A=A*V'; LM=LM/norm(LM(:));
A=reshape(A,[dp,d1,d2,d3,d4]);
v4=reshape( reshape(T,[dp*D^4,dp]) * reshape(A,[dp,d1*d2*d3*d4]),[dp,D,D,D,D,d1,d2,d3,d4] );
R=V*LM*V'; LM=diag(LM);

v4=reshape( permute(v4,[3,7,1,2,6,5,9,4,8]),[D*d2,dp*D*d1*D*d4*D*d3] );
v4=reshape( reshape( reshape(v2,[d2,D*d2])*v4,[d2*dp*D*d1*D*d4,D*d3] ) * (reshape(v3,[d3,D*d3]).'),[d2,dp,D*d1,D*d4,d3] );
v4=reshape( reshape(v1,[d1,D*d1]) * reshape( permute(v4,[3,2,1,5,4]),[D*d1,dp*d2*d3*D*d4] ),[d1*dp*d2*d3,D*d4] );
v4=reshape( reshape( permute(conj(A),[5,2,1,3,4]),[d4,d1*dp*d2*d3] )*v4,[d4,D,d4] );
Norm=norm(v4(:));
v4=v4/Norm;
end