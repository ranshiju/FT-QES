function [A1,lambda1,dc2,U,V]=RedefineA(A,d,dc,dc1,L,lambda,ISlambda)
Lh=round(L/2); lambda=diag(sqrt(lambda));

A=reshape( lambda^ISlambda*reshape( permute( reshape(A,[d^L,dc1,dc1]),[2,1,3] ),[dc1,d^L*dc1] ),[dc1*d^L,dc1] )*lambda^ISlambda;
A1=reshape( permute( reshape(A,[dc1,d^Lh,d^Lh,dc1]),[1,3,4,2] ),[dc1*d^Lh,dc1*d^Lh] );
[U,lambda1,V]=svds(A1,dc); S=size(lambda1); dc2=S(1);
U=reshape(U*sqrt(lambda1),[dc1,d^Lh*dc2]);
V=reshape(V*sqrt(lambda1),[dc1,d^Lh*dc2]);

A1=reshape( permute( reshape( V'*pinv(lambda^(ISlambda*2))*U,[d^Lh,dc2,d^Lh,dc2] ),[1,3,2,4] ),[d^L*dc2*dc2,1] );
lambda1=diag(lambda1); lambda1=lambda1/norm(lambda1);
end