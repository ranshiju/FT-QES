function [Ih,PatchH]=Finite2DHCcouplings(N)
if(N==16)
    Ih=[1,3;
            2,3;
            3,4;
            5,6;
            4,6;
            4,7;
            7,8;
            6,9;
            7,11;
            9,12;
            9,10;
            10,11;
            11,14;
            10,13;
            13,15;
            13,16];
        PatchH.h1=[]; % Patch h/3
        PatchH.h2=[1,2,5,8,12,14,15,16]; % Patch 2*h/3
elseif(N==8)
    Ih=[1,2;
            2,3;
            3,4;
            2,7;
            4,5;
            6,7;
            5,6;
            6,8];
        PatchH.h1=[1,8]; % Patch h/3
        PatchH.h2=[2,4,5,7]; % Patch 2*h/3
end
end