function [S,Ent,Xi,Err]=Entanglement_infinite_MPS(MPS,time,dtime,Eps,varargin)
% varargin{1}: ways of computing entanglement
% Way=1: QR decomposition
% Way=2: Eigs of the transfer matrix

time0=30;
Nv=numel(varargin);
if(Nv==0)
    Way=1;
elseif(Nv==1)
    Way=varargin{1};
end

Ni=numel(MPS);
if(Way==1)
    if(Ni==4)
    %     A1=MPS{1}; A2=MPS{2}; lamL=MPS{3}; lamR=MPS{4};
    %     A1=MPS_multiply_lambda_A(A1,sqrt(lamL),sqrt(lamR));
    %     A2=MPS_multiply_lambda_A(A2,sqrt(lamL),sqrt(lamR));
    elseif(Ni==1)
        A=MPS{1};
        [d,d1,d2]=size(A);
        A=permute(A,[2,1,3]); AL=A;
        AR=permute(AL,[3,2,1]);
        AL0=AL; AR0=AR;
    %     AR=AL;
        S0=ones(d1,1);

        for t=1:time
            [~,RL]=qr(reshape(AL,[d1*d,d2]),0);
            AL=RL*reshape(AL0,[d1,d*d2]);
            AL=AL/max(abs(AL(:)));

            [~,RR]=qr(reshape(AR,[d2*d,d1]),0);
            AR=RR*reshape(AR0,[d2,d*d1]);
            AR=AR/max(abs(AR(:)));

            if(rem(t,dtime)==0 && t>time0)
                [S]=svd(RL*(RR).');
                S=S/norm(S);
                Err=norm(S0-S);
                if(Err<Eps)
                    break
                else
                    S0=S;
                end
            end
        end
    end
    Xi=0;
elseif(Way==2)
    if(Ni==1)
        A=MPS{1};
        [d,d1,d2]=size(A);
        A=reshape(A,[d,d1*d2]);
        A=reshape( permute( reshape( A'*A,[d1,d2,d1,d2] ),[1,3,2,4] ),[d1*d1,d2*d2] );
        
        opt.maxit=10000;
        opt.v0=reshape( eye(d1),[d1*d1,1] );
        [vL,DL]=eigs(A,2,'LM',opt); vL=vL(:,1);
        
        opt.v0=reshape( eye(d2),[d2*d2,1] );
        [vR,DR]=eigs(A.',2,'LM',opt); vR=vR(:,1); 
        
        DL=abs(real(DL(1,1))/real(DL(2,2))); Xi(1)=1/log(DL);
        DR=abs(real(DR(1,1))/real(DR(2,2))); Xi(2)=1/log(DR);
        Xi=sum(Xi)/2;
        
        [~,~,S,~,~]=TruncationFromRM(reshape(vL,[d1,d1]),reshape(vR,[d2,d2]),ones(d1,1),d1);
        Err=0;
    end
end
Ent=Entanglement(S);
end