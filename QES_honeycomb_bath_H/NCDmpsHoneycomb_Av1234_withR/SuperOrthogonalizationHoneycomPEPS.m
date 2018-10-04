function [T,LM]=SuperOrthogonalizationHoneycomPEPS(T,LM,varargin)
Ni=numel(varargin);
if(Ni==1)
    dc=varargin{1}*ones(3,1);
else
    dc=zeros(1,3);
    for x=1:3
        dc(x)=numel(LM{x});
    end
end

Time=1000;
Time0=2;
ConvEps=1e-10;

for t=1:Time
    LM0=LM;
    for x=1:3
        LM1=LM;
        LM1{x}=ones(size(LM1{x}));
        ML=PEPS_bondRM(T{1},x,LM1);
        MR=PEPS_bondRM(T{2},x,LM1);
        [UL,UR,LM{x}]=TruncationFromRM(ML,MR,LM{x},dc(x));
        T{1}=TransforPEPS(T{1},UL,x);
        T{2}=TransforPEPS(T{2},UR,x);
    end
    
    if(t>min(Time0,2))
        Conv= norm(LM0{1}-LM{1})+norm(LM0{2}-LM{2})+norm(LM0{3}-LM{3});
        if(Conv<3*ConvEps)
            break;
        end
    end
end
end

function T=TransforPEPS(T,U,x)
S=size(T); Nb=numel(S); dd=size(U,2);
T=reshape( permute(T,[1:x,x+2:Nb,x+1]),[prod(S([1:x,x+2:Nb])),S(x+1)] ) * U;
T=permute( reshape(T,[S([1:x,x+2:Nb]),dd]),[1:x,Nb,x+1:Nb-1] );
end