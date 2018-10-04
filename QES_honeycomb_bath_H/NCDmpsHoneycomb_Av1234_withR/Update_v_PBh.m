function [vnew,An,R,LM,Norm]=Update_v_PBh(A,Tcore,U,Ht,Hb,v3,v4,bond,ISLM)
d=2; 

A=permute(A,[1:bond,bond+2:5,bond+1]);
[dp,d1,d2,d3,d4]=size(A);

if(ISLM)
    [A,LM,V]=svd(reshape(A,[dp*d1*d2*d3,d4]),0); 
    A=A*V'; LM=LM/norm(LM(:));
    R=V*LM*V'; LM=diag(LM);
else
    [A,R]=qr(reshape(A,[dp*d1*d2*d3,d4]),0); 
    LM=[];
end

A=reshape(A,[dp,d1,d2,d3,d4]);
An=permute(A,[1:bond,5,bond+1:4]);

[~,d1,d2,d3,d4]=size(An);
A=reshape(An,[d,d,d,d1,d2,d3,d4]);
Ac=Update_A_PBhamiltonian(conj(A),conj(Ht),[1,2,3]);

if(bond==3 || bond==4)
    A1=Update_A_PBhamiltonian(A,Hb{1},[4,1]);
    A1=Update_A_PBhamiltonian(A1,Hb{2},[5,2]);
    if(bond==3)
        Hbt_3=Hbt_half(Tcore,U,v3,v4,3);
        vnew=ncon({Ac,A1,Hbt_3},{[3,4,6,1,2,-1,5],[3,4,8,1,2,-3,7],[5,6,7,8,-2]});
    elseif(bond==4)
        Hbt_4=Hbt_half(Tcore,U,v3,v4,4);
        vnew=ncon({Ac,A1,Hbt_4},{[3,4,5,1,2,6,-1],[3,4,8,1,2,7,-3],[6,5,7,8,-2]});
    end
elseif(bond==1 || bond==2)
    U=reshape(U,[d,d,d*d]);
    A1=Update_A_PBhamiltonian(A,Hb{3},[6,7,3]);
    if(bond==1)
        A1=Update_A_PBhamiltonian(A1,Hb{2},[5,2]);
        vnew=ncon({Ac,A1,U},{[6,4,5,-1,1,2,3],[7,4,5,-3,1,2,3],[6,7,-2]});
    elseif(bond==2)
        A1=Update_A_PBhamiltonian(A1,Hb{1},[4,1]);
        vnew=ncon({Ac,A1,U},{[4,6,5,1,-1,2,3],[4,7,5,1,-3,2,3],[6,7,-2]});
    end
end

Norm=norm(vnew(:));
vnew=vnew/Norm;
end

function Hbt_34=Hbt_half(Tcore,U,v3,v4,bond)
[d1,D,~]=size(v3); d=2;

if(bond==3)
    v4=reshape( permute(v4,[2,1,3]),[D,d1*d1] );
    Hbt_34=Tucker_Product(Tcore,{eye(D),v4,U.'});
    Hbt_34=permute( reshape(Hbt_34,[D,d1,d1,d,d]),[2,4,3,5,1] );
elseif(bond==4)
    v3=reshape( permute(v3,[2,1,3]),[D,d1*d1] );
    Hbt_34=Tucker_Product(Tcore,{v3,eye(D),U.'});
    Hbt_34=permute( reshape(Hbt_34,[d1,d1,D,d,d]),[1,4,2,5,3] );
end
end