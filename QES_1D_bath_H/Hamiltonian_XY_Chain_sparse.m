function H=Hamiltonian_XY_Chain_sparse(Jx,hz,L)
d=2;

H0=sparse(XYmodel_Hamiltonian(Jx,hz));
I=sparse(eye(d,d));

for p=1:L-1
    H1=H0;
    for n=1:p-1
        H1=kron(H1,I);
    end
    for n=1:L-p-1
        H1=kron(I,H1);
    end
    
    if(p==1)
        H=H1;
    else
        H=H+H1;
    end
end
end