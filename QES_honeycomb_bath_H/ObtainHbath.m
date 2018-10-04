function [Hb]=ObtainHbath(HL,HR,d1,d2,D)
% Bath-physical Hamiltonian
HL=reshape( permute( reshape(HL,[d1,D,d1]),[1,3,2] ),[d1*d1,D] );
HR=reshape( permute( reshape(HR,[d2,D,d2]),[2,1,3] ),[D,d2*d2] );
Hb=HL*HR;
Hb=reshape( permute( reshape(Hb,[d1,d1,d2,d2]),[1,3,2,4] ),[d1*d2,d1*d2] );
end