function [Hbt1,Hbt2,dg]=Boundary_Interactions_triangle(U,Tcore,v3,v4,ISsvd)
[d1,D,~]=size(v3); d=2;

v3=reshape( permute(v3,[2,1,3]),[D,d1*d1] );
v4=reshape( permute(v4,[2,1,3]),[D,d1*d1] );

if(ISsvd)
    Eps=1e-15;
    
    Hbt=Tucker_Product(Tcore,{v3,v4,U.'});
    Hbt=permute( reshape(Hbt,[d1,d1,d1,d1,d,d]),[1,3,2,4,5,6] );

    [Hbt1,tmp,Hbt2]=svd( reshape(Hbt,[d1^4,d^2]),0 );
    dg=find(diag(tmp)>Eps,1,'last');
    Hbt1=reshape( Hbt1(:,1:dg)*sqrt(tmp(1:dg,1:dg)),[d1*d1,d1*d1,dg] );
    Hbt2=reshape( conj(Hbt2(:,1:dg))*sqrt(tmp(1:dg,1:dg)),[d,d,dg] );
else
    Hbt1=Tucker_Product(Tcore,{v3,v4,eye(size(U.'))});
    Hbt1=reshape( permute( reshape(Hbt1,[d1,d1,d1,d1,d*d]),[1,3,2,4,5] ),[d1*d1,d1*d1,d*d] );
    Hbt2=reshape( U,[d,d,d*d] );
    dg=d*d;
end
end

