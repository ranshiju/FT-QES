function [v,lambda,X]=Orthogonalize_v(v,D,dc)
[U,lambda,V]=svd(reshape(v,[dc,D*dc]),0);
V=V(:,1:dc); lambda=diag(lambda);
v=reshape( U*V',[dc*D*dc,1] );
X=U*diag(lambda(1:dc))*U';
end