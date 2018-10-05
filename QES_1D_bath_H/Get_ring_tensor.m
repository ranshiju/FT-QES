function M=Get_ring_tensor(Au,Ad,vL,vR,d,D,dc)
M=reshape( permute( reshape( reshape(permute(vL,[1,3,2]),[D*dc,dc])*reshape( permute(Au,[2,1,3]),[dc,d*dc] ),[D,dc,d,dc] ),[3,1,4,2] ),[d*D,dc*dc] );
M=M*reshape( permute( reshape( reshape(Ad,[d*dc,dc])*reshape(vR,[D*dc,dc]).',[d,dc,D,dc] ),[4,2,1,3] ),[dc*dc,d*D] );
M=reshape(M,[d,D,d,D]);
end