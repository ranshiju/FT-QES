function H=FullHamiltonian(Hbath,Hbulk,L)
d=2;
D=size(Hbath{1},1); D=round(D/d);

% O2=sparse(zeros(d,d));
% H0=O2;
% for n=1:L-2
%     H0=kron(H0,O2);
% end
% H=kron(H0,O2);
H=sptensor([],[],[D,ones(1,L-2)*d,D,D,ones(1,L-2)*d,D]);

for n=1:L-1
    if(n==1)
        H1=Hbath{1}; d1=D; d2=d;
        Did=d^(L-3)*D; Sid=[ones(1,L-3)*d,D,ones(1,L-3)*d,D];
    elseif(n==L-1)
        H1=Hbath{2}; d1=d; d2=D;
        Did=d^(L-3)*D; Sid=[D,ones(1,L-3)*d,D,ones(1,L-3)*d];
    else
        H1=Hbulk; d1=d; d2=d;
        Did=d^(L-4)*D^2; Sid=[D,ones(1,L-4)*d,D,D,ones(1,L-4)*d,D];
    end
    Hid=reshape( sptensor([(1:Did)',(1:Did)'],ones(Did,1),[Did,Did]),Sid );
    H1=reshape( sptensor(H1),[d1,d2,d1,d2] );
    
    % Pertmute Hbath to the right position
    V=[1,2,5:L+2,3,4,L+3:2*L];
    V=ReOrderV(V,[3:n+1,1,2,n+2:L,[3:n+1,1,2,n+2:L]+L]);

    H=H+permute(ttt(H1,Hid),V);
end
H=reshape(H,[D*d^(L-2)*D,D*d^(L-2)*D]);
end