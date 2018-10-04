function Ort=TestSuperOrthogonality(T,LM)

Ort=0;
for x=1:3
    ML=PEPS_bondRM(T{1},x,LM);
    MR=PEPS_bondRM(T{2},x,LM);
    
    ML=diag(LM{x}.^(-1))*ML*diag(LM{x}.^(-1));
    ML=ML/ML(1);
    ML=eye(size(ML))-ML;
    Ort=Ort+norm(ML(:))/6;
    
    MR=diag(LM{x}.^(-1))*MR*diag(LM{x}.^(-1));
    MR=MR/MR(1);
    MR=eye(size(MR))-MR;
    Ort=Ort+norm(MR(:))/6;
end
end