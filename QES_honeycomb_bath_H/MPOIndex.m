function [Index,CoupPhys,PosP,PosB]=MPOIndex(Cluster)
% Positions MUST be in ascending order
if Cluster==1
    Index.couplings{1}=[4,5];
    Index.positions{1}=[1,2,3];

    Index.couplings{2}=[3]; %#ok<*NBRAK>
    Index.positions{2}=[2,6];

    Index.couplings{3}=[4,2,1,3,2,1,7];
    Index.positions{3}=[4,5,6,7,9,11,12,13];

    Index.couplings{4}=[5];
    Index.positions{4}=[7,8];

    Index.couplings{5}=[6];
    Index.positions{5}=[9,10];

    Index.couplings{6}=[3];
    Index.positions{6}=[5,12];

    Index.couplings{7}=[3,7];
    Index.positions{7}=[11,14,15];

    Index.couplings{8}=[6];
    Index.positions{8}=[14,16];

    % Physical couplings (three kinds), for obvervation
    CoupPhys=... 
           [2,6,3;
            5,6,2;
            6,7,1;
            5,12,3;
            7,9,3;
            11,12,1;
            9,11,2;
            11,14,3];
    
    Index.MidPhysBond = 1:8;
    PosP=[2,5,6,7,9,11,12,14]; % Positions of physical sites
    PosB=[1,3,4,8,10,13,15,16]; % Positions of bath sites
        
elseif Cluster==2
    n = 0;
    n = n+1;
    Index.positions{n}=[1,2,3];
    Index.couplings{n}=[4,5];
    
    n = n+1;
    Index.positions{n}=[4,5,6];
    Index.couplings{n}=[4,5];
    
    n = n+1;
    Index.positions{n}=[2,9];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[5,11];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[7:12,14,16:20];
    Index.couplings{n}=[4,2,1,2,1,3,2,1,2,1,7];
    
    n = n+1;
    Index.positions{n}=[12,13];
    Index.couplings{n}=[5];
    
    n = n+1;
    Index.positions{n}=[14,15];
    Index.couplings{n}=[6];
    
    n = n+1;
    Index.positions{n}=[10,17];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[8,19];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[18,21,22];
    Index.couplings{n}=[3,7];
    
    n = n+1;
    Index.positions{n}=[16,24,25];
    Index.couplings{n}=[3,6];
    
    n = n+1;
    Index.positions{n}=[21,23];
    Index.couplings{n}=[6];
    
    n = n+1;
    Index.positions{n}=[24,26];
    Index.couplings{n}=[6];
    
    CoupPhys=... 
           [2,9,3;
            5,11,3;
            8,9,2;
            9,10,1;
            10,11,2;
            11,12,1;
            12,14,3;
            8,19,3;
            10,17,3;
            14,16,2;
            16,17,1;
            17,18,2;
            18,19,1;
            16,24,3;
            18,21,3
            ];
        Index.MidPhysBond = [4, 5, 9, 11, 12];
        
        PosP=[2,5,8,9,10,11,12,14,16,17,18,19,21,24];
        PosB=[1,3,4,6,7,13,15,20,22,23,25,26];
        
elseif Cluster==3
    n = 1;
    Index.positions{n}=[1,2,3];
    Index.couplings{n}=[4,5];
    
    n = n+1;
    Index.positions{n}=[2, 6, 7, 8];
    Index.couplings{n}=[3, 1, 5];
    
    n = n+1;
    Index.positions{n}=[4, 5, 6];
    Index.couplings{n}=[4, 2];
    
    n = n+1;
    Index.positions{n}=[5, 13];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[7, 11];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[9 ,11, 12, 13, 15, 16, 18, 19, 20, 21, 22];
    Index.couplings{n}=[1, 2, 1, 2, 3, 1, 2, 1, 2, 6];
    
    n = n+1;
    Index.positions{n}=[9, 10];
    Index.couplings{n}=[5];
    
    n = n+1;
    Index.positions{n}=[9, 21];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[12, 19];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[14, 15];
    Index.couplings{n}=[4];
    
    n = n+1;
    Index.positions{n}=[16, 17];
    Index.couplings{n}=[7];
    
    n = n+1;
    Index.positions{n}=[18, 26];
    Index.couplings{n}=[3];
    
    n = n+1;
    Index.positions{n}=[20, 23, 24];
    Index.couplings{n}=[3, 6];
    
    n = n+1;
    Index.positions{n}=[23, 25, 26, 27];
    Index.couplings{n}=[2, 1, 7];
    
    n = n+1;
    Index.positions{n}=[25, 28, 29];
    Index.couplings{n}=[3, 7];
    
    n = n+1;
    Index.positions{n}=[28, 30];
    Index.couplings{n}=[6];
    
    CoupPhys=... 
           [2,6,3;
            5,6,2;
            6,7,1;
            5,13,3;
            7,11,3;
            11, 12, 2;
            12, 13, 1;
            9, 11, 1;
            13, 15, 2;
            15, 16, 3;
            12, 19, 3;
            9, 21, 3;
            20, 21, 2;
            19, 20, 1;
            18, 19, 2;
            16, 18, 1;
            18, 26, 3;
            20, 23, 3;
            23, 25, 2;
            25, 26, 1;
            25, 28, 3];
        Index.MidPhysBond = [6, 7, 11, 14, 15];
        
    PosB=[1, 3, 4, 8, 10, 14, 17, 22, 24, 27, 29, 30];
    PosP = GetPosPfromPosB(30, PosB);
end

Index.Ntwobody = 0;
for n = 1:numel(Index.couplings)
    Index.Ntwobody = Index.Ntwobody + numel(Index.couplings{n});
end
end

function PosP = GetPosPfromPosB(Ntot, PosB)
PosP = zeros(1, Ntot-numel(PosB));
kk = 0;
for k = 1:Ntot
    if ~ismember(k, PosB)
        kk = kk +1;
        PosP(kk) = k;
    end
end
end