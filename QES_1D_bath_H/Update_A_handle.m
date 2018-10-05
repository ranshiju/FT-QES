function [A1]=Update_A_handle(A,v,G,para)
d=2;
if(para.Otrotter==2)
    A=Evolve_H_Bulk(A,G.bulk,para.L,2,para.dc1,para.ISsparse);
end
[A1]=Update_A_BoundaryH(A,v,G.L,G.R,para.L,para.dc1);
A1=reshape( Evolve_H_Bulk(A1,G.bulk,para.L,2,para.dc1,para.ISsparse),[d^para.L*para.dc*para.dc,1] );
end