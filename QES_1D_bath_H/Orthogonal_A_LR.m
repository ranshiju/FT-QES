function [AL,AR,lambdaL,lambdaR]=Orthogonal_A_LR(A,d,dc,L)
AL=reshape(A,[d^L*dc,dc]);
[AL,lambdaL,~]=svd(AL,0);
% [AL,~]=qr(AL,0);
AL=reshape(AL,[d^L,dc,dc]);

AR=reshape( permute(reshape(A,[d^L,dc,dc]),[1,3,2]),[d^L*dc,dc] );
[AR,lambdaR,~]=svd(AR,0);
% [AR,~]=qr(AR,0);
AR=permute( reshape(AR,[d^L,dc,dc]),[1,3,2] );

lambdaL=diag(lambdaL); lambdaR=diag(lambdaR);
end