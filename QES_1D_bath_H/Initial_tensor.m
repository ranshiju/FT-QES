function [A,v]=Initial_tensor(L,dc,d,Way)
% v.L and v.R: (dc,d*d,dc)

A1=randn([ones(1,L)*d,dc,dc]);
A=zeros(size(A1));
for n=1:L
    A1=permute(A1,[2:L,1,L+1,L+2]);
    A=A+A1;
end
A=reshape( A+permute(A,[1:L,L+2,L+1]),[d^L*dc*dc,1] );
A=A/norm(A(:));

if(Way==1)
    D=d;
elseif(Way==2)
    D=d*d;
end
% v.L=reshape( permute( reshape( eye(dc*D),[dc,D,dc,D] ),[1,2,4,3] ),[dc*D*D*dc,1] );
v.L=rand(dc*D,dc*D);
v.L=reshape( permute( reshape(v.L+v.L',[dc,D,dc,D]),[2,1,3,4] ),[dc*D*D*dc,1] );

v.L=v.L/norm(v.L); v.R=v.L;
end