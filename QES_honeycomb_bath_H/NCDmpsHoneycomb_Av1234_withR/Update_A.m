function A1=Update_A(T,A,v1,v2,v3,v4,S)
[dp,D,~,~,~,~]=size(T);
d1=S(1); d2=S(2); d3=S(3); d4=S(4);

A1=reshape( reshape(T,[dp*D^4,dp]) * reshape(A,[dp,d1*d2*d3*d4]),[dp,D,D,D,D,d1,d2,d3,d4] );
A1=reshape( permute(A1,[3,7,1,2,6,5,9,4,8]),[D*d2,dp*D*d1*D*d4*D*d3] );
A1=reshape( reshape( reshape(v2,[d2,D*d2])*A1,[d2*dp*D*d1*D*d4,D*d3] ) * (reshape(v3,[d3,D*d3]).'),[d2,dp,D*d1,D*d4,d3] );
A1=reshape( reshape(v1,[d1,D*d1]) * reshape( permute(A1,[3,2,1,5,4]),[D*d1,dp*d2*d3*D*d4] ),[d1*dp*d2*d3,D*d4] );
A1=reshape( permute( reshape( A1 * (reshape(v4,[d4,D*d4]).'),[d1,dp,d2,d3,d4] ),[2,1,3,4,5] ),[dp*d1*d2*d3*d4,1] );
end