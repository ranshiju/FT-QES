function [UL,UR,lambda,dc,Norm]=TruncationFromRM(ML,MR,lambda,dc,varargin)
% [UL,UR,lambda,dc,Norm]=TruncationFromRM(ML,MR,lambda,varargin)
% Get truncation matrices from reduced matrices
% varargin{1} judges if lambda would be kept
% varargin{2} give Eps to truncate very small eigenvalues

if(numel(varargin)==1)
    ExistLambda=varargin{1}; Eps=0; ISEIG=0;
elseif(numel(varargin)==2)
    ExistLambda=varargin{1}; Eps=varargin{2}; ISEIG=0;
elseif(numel(varargin)==3)
    ExistLambda=varargin{1}; Eps=varargin{2}; ISEIG=varargin{3};
else
    ExistLambda=1; Eps=0; ISEIG=0;
end

% [UL,DL,~]=svd(ML);
[UL,DL]=eig(ML,'nobalance'); %DL=abs(DL);
% [UL,DL]=Inverse_EigV(UL,diag(DL));
[UL,DL]=Eig_In_Order_Abs(UL,DL); DL=diag(DL);
DL=DL/norm(DL);

% [UR,DR,~]=svd(MR);
[UR,DR]=eig(MR,'nobalance'); %DR=abs(DR); 
% [UR,DR]=Inverse_EigV(UR,diag(DR));
[UR,DR]=Eig_In_Order_Abs(UR,DR); DR=diag(DR);
DR=DR/norm(DR);

% dc=max(size(DL));
dc1=min( [find(abs(DL)>Eps,1,'last'), find(abs(DR)>Eps,1,'last')]);
if(dc1<1.1)
    warning('dc<=1 in TruncationRM. \n')
    if(DL(2)<5*Eps)
        pinvDL=[1,0];
    else
        pinvDL=sqrt(DL(1:2)).^(-1);
    end
    if(DR(2)<2*Eps)
        pinvDR=[1,0];
    else
        pinvDR=sqrt(DR(1:2)).^(-1);
    end
    dc1=2;
    warning('Only one large D. \n')
else
    pinvDL=sqrt(DL(1:dc1)).^(-1); 
    pinvDR=sqrt(DR(1:dc1)).^(-1);
end

if(numel(DL)~=dc1)
%     fprintf('Some eigenvalues are too small with d1 = %g, d2 = %g \n',S,dc)
    UL=UL(:,1:dc1); UR=UR(:,1:dc1);
    DL=DL(1:dc1); DR=DR(1:dc1);
end

DL=sqrt(DL); DR=sqrt(DR);

if(ISEIG)
    [U,lambda]=eig( diag(DL)*UL'*diag(lambda)*conj(UR)*diag(DR) );
    [U,lambda]=Eig_In_Order_Abs(U,lambda);
    V=conj(U);
else
    [U,lambda,V]=svd( diag(DL)*UL'*diag(lambda)*conj(UR)*diag(DR) );
end

dc=min(dc,dc1);
lambda=diag(lambda(1:dc,1:dc)); Norm=norm(lambda);
lambda=lambda/Norm;
% dc=min( find(lambda>Eps,1,'last'),dc );
% U=U(:,1:dc); V=V(:,1:dc); lambda=lambda(1:dc);

if(ExistLambda)
    UL=UL*diag(pinvDL)*U(:,1:dc)*sqrt(Norm);
    UR=UR*conj(diag(pinvDR))*conj(V(:,1:dc))*sqrt(Norm);
else
    UL=UL*diag(pinvDL)*U(:,1:dc)*diag(sqrt(lambda))*sqrt(Norm);
    UR=UR*diag(pinvDR)*conj(V(:,1:dc))*diag(sqrt(lambda))*sqrt(Norm);
end
Norm=1;
end