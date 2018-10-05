function A=Evolve_H_Bulk(A,Gbulk,L,d,dc,ISsparse)
% Input: (d,d,...,d,dc,dc)
% Output: (d,d,...,d,dc,dc)

if(ISsparse)
    A=Gbulk*reshape(A,[d^L,dc*dc]);
else
%     Time=round(L/2);
%     for t=1:Time
%         A=reshape( G.G2*reshape(A,[d*d,d^(L-2)*dc*dc]),[ones(1,L)*d,dc,dc] );
%         A=permute(A,[3:L,1,2,L+1,L+2]);
%     end
% 
%     A=permute(A,[2:L,1,L+1,L+2]);
%     for t=1:Time-1
%         A=reshape( G.G2*reshape(A,[d*d,d^(L-2)*dc*dc]),[ones(1,L)*d,dc,dc] );
%         A=permute(A,[3:L,1,2,L+1,L+2]);
%     end
end
% if(Way==1)
%     A=permute(A,[L,L+1,1:L-2,L+2,L-1]); 
% end
end