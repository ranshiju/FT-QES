function [A,lambda,dc2,Err,Canon]=Canonicalize_A(A,d,dc1,dc,time)
% [UL,UR]=TM_MPS_eig(A,d,dc,L);
% [vL,vR,~,dc1]=TruncationFromRM(UL,UR,ones(dc,1),0);
% 
% I=vL*vR.'; I=I/I(1); I=I-eye(dc1);
% Err=norm(I(:))/dc1;
% 
% A1=reshape( reshape( vL.'*reshape( permute( reshape(A,[d^L,dc,dc]),[2,1,3] ),[dc,d^L*dc] ),[dc1*d^L,dc] )*vR,[dc1,d^L,dc1] );
% A1=reshape( permute(A1,[2,1,3]),[d^L*dc1*dc1,1] );

Eps=1e-12; Canon=1;
lambda0=ones(dc1,1); lambda=lambda0;
vL=eye(dc1); vR=vL; dd=dc1; I=eye(dc1);
A=permute(A,[2,1,3]); 
for t=1:time
    UL=reshape( diag(lambda)*vR.'*reshape( permute(A,[3,2,1]),[dc1,d*dc1] ),[dd*d,dc1] )*vL;
    UL=UL'*UL;
    UR=reshape( diag(lambda)*vL.'*reshape( A,[dc1,d*dc1] ),[dd*d,dc1] )*vR;
    UR=UR'*UR;

    [UL,UR,lambda,dd]=TruncationFromRM(UL,UR,lambda,1,1e-15);
    vL=vL*UL; vR=vR*UR;
    vL=vL/max(abs(vL(:))); vR=vR/max(abs(vR(:)));
%     A=reshape( reshape( vL.'*reshape( A,[dc,d^L*dc] ),...
%         [dc*d^L,dc] )*vR,[dc,d^L,dc] );

    UL=abs(UL/UL(1))-I; UR=abs(UR/UR(1))-I;
    Canon=(norm(UL(:))+norm(UR(:)))/2/dc1;

    if(numel(lambda)==numel(lambda0))
        Err=norm(lambda-lambda0);
    else
        Err=1;
    end
%     vL=vL/vL(1); vR=vR/vR(1);
%     Err=abs(vL-eye(dc))+abs(vR-eye(dc));
%     Err=norm(Err(:))/2/dc;
    if(Err<Eps && Canon<1e-8)
%         keyboard
        break;
    else
        lambda0=lambda;
    end
%     keyboard
end
if(Err>Eps || Canon>1e-8)
    warning('Bad convergence in canonicalization: Er=%g  ErCan=%g \n',Err,Canon);
end

dc2=min(dd,dc);
vL=vL(:,1:dc2)*diag(sqrt(lambda(1:dc2))); vR=vR(:,1:dc2)*diag(sqrt(lambda(1:dc2)));
A=permute( reshape( reshape( vL.'*reshape(A,[dc1,d*dc1]),[dc2*d,dc1] )...
    *vR,[dc2,d,dc2] ),[2,1,3] );

% vL=permute( reshape( vL,[d*d,round(dc1/d/d),dc2] ),[3,1,2] ); 
% vR=permute( reshape( vR,[d*d,round(dc1/d/d),dc2] ),[3,1,2] ); 
% test code ------------------------------------------------
% A1=reshape(A,[d^L,dc*dc]);
% M=reshape( permute( reshape( A1'*A1,[dc,dc,dc,dc] ),[1,3,2,4] ),[dc*dc,dc*dc] );
% v=reshape( diag(lambda),[dc*dc,1] ); v=v/norm(v);
% v1=v.'*M; v1=v1.'/norm(v1);
% ErrEig(1)=norm(v1-v);
% v1=M*v; v1=v1/norm(v1);
% ErrEig(2)=norm(v1-v);
% keyboard
end