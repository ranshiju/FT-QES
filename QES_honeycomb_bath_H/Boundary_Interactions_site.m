function [Hb1,Hb2,dg, Hb]=Boundary_Interactions_site(v,U,ISsvd)
[d1,D,~]=size(v); d=2;

Hb=reshape( permute(v,[1,3,2]),[d1*d1,D] ) * U.';
if(ISsvd)
    Eps=1e-15;
    
    [Hb1,tmp,Hb2]=svd(Hb,0);
    dg=find(diag(tmp)>Eps,1,'last');
    Hb1=reshape( Hb1(:,1:dg)*sqrt(tmp(1:dg,1:dg)),[d1,d1,dg] );
    Hb2=reshape( conj(Hb2(:,1:dg))*sqrt(tmp(1:dg,1:dg)),[d,d,dg] );
else
    Hb1=permute(v,[1,3,2]);
    Hb2=reshape(U,[d,d,d*d]);
    dg=d*d;
end
Hb1=permute(Hb1,[1,3,2]); Hb2=permute(Hb2,[1,3,2]);
Hb = Hb - trace(Hb)*eye(size(Hb))/d1/d;
end