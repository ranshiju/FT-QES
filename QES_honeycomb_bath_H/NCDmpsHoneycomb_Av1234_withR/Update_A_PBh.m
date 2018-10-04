function A1=Update_A_PBh(A,Ht,Hb,dc)
d=2; d1=dc; d2=dc; d3=dc; d4=dc;
% [dp,d1,d2,d3,d4]=size(A);

A1=reshape(A,[d,d,d,d1,d2,d3,d4]);
A1=Update_A_PBhamiltonian(A1,Hb{1},[4,1]);
A1=Update_A_PBhamiltonian(A1,Hb{2},[5,2]);
A1=Update_A_PBhamiltonian(A1,Hb{3},[6,7,3]);
A1=Update_A_PBhamiltonian(A1,Ht,[1,2,3]);
A1=reshape(A1,[d^3*d1*d2*d3*d4,1]);
end