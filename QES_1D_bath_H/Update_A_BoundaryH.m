function [A1]=Update_A_BoundaryH(A,v,GL,GR,L,dc)
% Input & Output: (d,d,...,d,dc,dc)
d=2; A=permute(reshape(A,[ones(1,L)*d,dc,dc]),[1,L+1,2:L-1,L+2,L]); % (d,dc,d,...,d,dc,d)

% if(Way==1)
GL=reshape( reshape( permute( reshape(v.R,[dc,d*d,dc]),[1,3,2] ),[dc*dc,d*d] ) * ...
    reshape( permute( reshape(GL,[d,d*d,d]),[2,1,3] ),[d*d,d*d] ),[dc,dc,d,d] );
GR=reshape( reshape( permute( reshape(v.L,[dc,d*d,dc]),[1,3,2] ),[dc*dc,d*d] ) * ...
    reshape( permute( reshape(GR,[d,d*d,d]),[2,1,3] ),[d*d,d*d] ),[dc,dc,d,d] );

GL=reshape( permute(GL,[2,4,1,3]),[dc*d,dc*d] );
GR=reshape( permute(GR,[3,1,4,2]),[d*dc,d*dc] );

A1=reshape( permute( reshape( reshape( GR*reshape(A,[d*dc,d^(L-2)*dc*d]),[d*dc*d^(L-2),dc*d] )*GL,[d,dc,d^(L-2),dc,d] ),[1,3,5,2,4] ),[d^L*dc*dc,1] );
% elseif(Way==2)
%     D=d^4;
% 
%     v.L=reshape( permute( reshape(v.L,[dc,D,dc]),[2,1,3] ),[D,dc*dc] );
%     v.R=reshape( permute( reshape(v.R,[dc,D,dc]),[2,1,3] ),[D,dc*dc] );
%     GL=reshape( permute( reshape( reshape( permute( reshape(GL*reshape(GL,[d,d*d*d]),[d,d*d,d*d,d]),[1,4,2,3] ),...
%         [d*d,d^4] )*v.L,[d,d,dc,dc] ),[4,2,3,1] ),[d*dc,d*dc] );
%     GR=reshape( permute( reshape( reshape( permute( reshape(GR*reshape(GR,[d,d*d*d]),[d,d*d,d*d,d]),[1,4,2,3] ),...
%         [d*d,d^4] )*v.R,[d,d,dc,dc] ),[1,3,2,4] ),[d*dc,d*dc] );
% 
%     A1=reshape( permute( reshape( reshape( GR*reshape(A,[d*dc,d^(L-2)*dc*d]),[d*dc*d^(L-2),dc*d] )*GL,[d,dc,d^(L-2),dc,d] ),[1,3,5,2,4] ),[d^L*dc*dc,1] );
% end
end