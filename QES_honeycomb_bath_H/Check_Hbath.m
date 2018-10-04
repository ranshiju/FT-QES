function H=Check_Hbath(v,U)
[chi,D,~]=size(v);
d=size(U,1);

v=reshape( permute(v,[1,3,2]),[chi*chi,D] );
U=reshape( permute(U,[1,3,2]),[d*d,D] );

H=reshape( permute( reshape( v*(U.'),[chi,chi,d,d] ),[1,3,2,4] ),[chi*d,chi*d] );
end