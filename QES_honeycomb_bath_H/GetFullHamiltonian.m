function H=GetFullHamiltonian(Ih2,Hp,Hb,PositionP,PositionB,Dim)
Np=size(Ih2.Hp,1); 
Nb=size(Ih2.Hb,1);

Lp=numel(PositionP);
Lb=numel(PositionB);
L=Lp+Lb;
d=Dim(PositionP(1));
D=Dim(PositionB(1));

H=sptensor([],[],[Dim,Dim]);
Hp=reshape( sptensor(Hp),[d,d,d,d] );

% Add physical interactions
Did=d^(Lp-2)*D^Lb;
for n=1:Np    
    DimId=Dim; DimId(Ih2.Hp(n,:))=[];
    Hid=reshape( sptensor([(1:Did)',(1:Did)'],ones(Did,1),[Did,Did]),[DimId,DimId] );
    
    V=[1,2,5:L+2,3,4,L+3:2*L];
    V=ReOrderV(V,[3:Ih2.Hp(n,1)+1,1,Ih2.Hp(n,1)+2:Ih2.Hp(n,2),2,Ih2.Hp(n,2)+1:L,[3:Ih2.Hp(n,1)+1,1,Ih2.Hp(n,1)+2:Ih2.Hp(n,2),2,Ih2.Hp(n,2)+1:L]+L]);
    
    H=H+reshape( permute( ttt(Hp,Hid),V ),[Dim,Dim] );
end

% Add physical-bath interactions
Did=d^(Lp-1)*D^(Lb-1);
for n=1:Nb
    if(Ih2.Hb(n,4)==0)
        Hnow=reshape(Hb{Ih2.Hb(n,3)},[d,D,d,D]);
    else
        Hnow=reshape(Hb{Ih2.Hb(n,3)},[D,d,D,d]);
    end
    Hnow=sptensor(Hnow);
    
    DimId=Dim; DimId(Ih2.Hb(n,1:2))=[];
    Hid=reshape( sptensor([(1:Did)',(1:Did)'],ones(Did,1),[Did,Did]),[DimId,DimId] );
    
    V=[1,2,5:L+2,3,4,L+3:2*L];
    V=ReOrderV(V,[3:Ih2.Hp(n,1)+1,1,Ih2.Hp(n,1)+2:Ih2.Hp(n,2),2,Ih2.Hp(n,2)+1:L,[3:Ih2.Hp(n,1)+1,1,Ih2.Hp(n,1)+2:Ih2.Hp(n,2),2,Ih2.Hp(n,2)+1:L]+L]);
    
    H=H+reshape( permute( ttt(Hnow,Hid),V ),[Dim,Dim] );
end
H=reshape(H,[d^Lp*D^Lb,d^Lp*D^Lb]);
end