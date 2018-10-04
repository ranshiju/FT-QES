function [EnvL,EnvR]=ObtainEnvsMPS(A)
% MPS has open boundary condition
% EnvL{1} and EnvR{end} are empty
%             2
%              |  
%  1   -    A     -  3

N=numel(A);
EnvL=cell(N,1); EnvR=cell(N,1); 

EnvL{2}=A{1}'*A{1};
for n=2:N-1
    EnvL{n+1}=Left2right(A{n},EnvL{n});
end

EnvR{N-1}=conj(A{N})*A{N}.';
for n=N-1:-1:2
    EnvR{n-1}=Right2left(A{n},EnvR{n});
end
end