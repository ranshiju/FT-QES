function v1=Right2left(A,v)
if(numel(v)==0)
    v1=conj(A)*(A.');
else
    [d1,d,d2]=size(A);
    v1=reshape( reshape(conj(A),[d1*d,d2])*v,[d1,d*d2] ) * (reshape(A,[d1,d*d2]).');
end
end