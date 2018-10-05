function H=Hamiltonian_Chain_sparse(H2,L)
d=2;

H0=sparse(H2);
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