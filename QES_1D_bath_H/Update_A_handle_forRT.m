function [A1]=Update_A_handle_forRT(A,v,G,Otrotter,L,dc,ISsparse)
d=2;
if(Otrotter==2)
    A=Evolve_H_Bulk(A,G.bulk,L,2,dc,ISsparse);
end
[A1]=Update_A_BoundaryH(A,v,G.L,G.R,L,dc);
A1=reshape( Evolve_H_Bulk(A1,G.bulk,L,2,dc,ISsparse),[d^L*dc*dc,1] );
end