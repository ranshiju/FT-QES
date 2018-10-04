function I0=Complement_vector(I1,Iall)
% 不得有重复数字
N=numel(I1);
I0=Iall;

for n=1:N
    x=find(I0==I1(n),1,'first');
    I0(x)=[];
end
end