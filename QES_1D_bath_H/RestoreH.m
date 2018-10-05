function H=RestoreH(J)
Spin={[0,0.5;0.5,0],[0,1i*0.5;-1i*0.5,0],[0.5,0;0,-0.5],eye(2)};
H=zeros(4,4);
for n1=1:4
    for n2=1:4
        H=H+J(n1,n2)*kron(Spin{n2},Spin{n1});
    end
end
end