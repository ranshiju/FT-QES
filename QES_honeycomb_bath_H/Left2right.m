function v1=Left2right(A,v)
if(numel(v)==0)
    v1=A'*A;
else
    [d1,d,d2]=size(A);
    v1=reshape(A,[d1*d,d2])' * reshape( v*reshape(A,[d1,d*d2]),[d1*d,d2] );
end
end