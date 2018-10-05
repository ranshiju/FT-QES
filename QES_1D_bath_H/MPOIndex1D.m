function [Index,CoupPhys, PosP, PosB]=MPOIndex1D(L)
% Positions MUST be in ascending order
% Index=1: physical (the 1st row of G)
% Index=n+1: n-th physical-bath or bath-physical (the (n+1)th row of G) [original order of physical and bath sites]
% Index=n+Nb+1: n-th bath-physical or physical-bath (the (n+5)th row of G) [the order is reversed]

Index.couplings{1}=[2,ones(1,L-3),3];
Index.positions{1}=1:L;

if rem(L, 2)==0
    Index.MidPhysBond = round(L/2):rand(L/2)+1;
else
    Index.MidPhysBond = round(L/2-0.2):round(L/2-0.2)+2;
end

PosP=2:L-1; % Positions of physical sites
PosB=[1, L]; % Positions of bath sites

CoupPhys=[(2:L-2).', (3:L-1).', ones(L-3, 1)];
end