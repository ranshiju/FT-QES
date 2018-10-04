function A1=Update_A_PBhamiltonian(A,H,Index)
[Ord,S]=Tensor_order(A);
Porder=Complement_vector(Index,1:Ord);
N0=prod(S(Index)); N1=prod(S(Porder));

Porder=[Index,Porder];
[~,Pback]=sort(Porder);
A1=permute(A,Porder);
S1=size(A1);
A1=permute( reshape( H*reshape(A1,[N0,N1] ),S1 ),Pback );
end