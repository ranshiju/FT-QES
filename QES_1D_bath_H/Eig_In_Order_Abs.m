function [U,D]=Eig_In_Order_Abs(U,D)
D1=D; U1=U; Dd=abs(diag(D));
[~,Order]=sort(Dd,'descend');
for m=1:numel(Dd)
    D(m,m)=D1(Order(m),Order(m));
    U(:,m)=U1(:,Order(m));
end

% p=max(size(D));
% UL1=zeros(p,p);DL1=zeros(1,p);
% SL=abs(diag(D));
% 
% for m1=p:-1:1
%     [~,Order]=max(SL);
%     DL1(1,p-m1+1)=D(Order,Order);
%     UL1(:,p-m1+1)=U(:,Order);
%     SL(Order)=-1;
% end
% 
% U=UL1; D=diag(DL1);
end