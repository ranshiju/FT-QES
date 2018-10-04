function expM=ExpmSptensor(M,tau,varargin)
% expm(tau*M) with M a sptensor and tau a small number (<1)
if(numel(varargin)==1)
    N=varargin{1};
else
    N=20; %Order of the Tylor expanssion
end

D=size(M,1);
expM=sptensor([(1:D).',(1:D).'],ones(D,1),[D,D]);
expM=expM+tau*M;

M1=M;
for n=2:N
    C=n*log(tau)-sum(log(1:n));
    M1=ttt(M1,M,2,1);
    expM=expM+exp(C)*M1;
end
end