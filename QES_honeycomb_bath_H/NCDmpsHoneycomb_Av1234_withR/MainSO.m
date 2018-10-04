function MainSO(varargin)
%% Preperation
fprintf('Start simulation of the physical-bath interactions \n');
d=2;
Ni=numel(varargin);
if(Ni==0)
    paraSO=ParameterSO;
    paraSO.PATHtree=[pwd,filesep];
    paraSO.EXP=['ODTNS_J_(',num2str(paraSO.Jh),',',num2str(paraSO.Jk),')_h_(',num2str(paraSO.hx),',',num2str(paraSO.hz),')_dc',num2str(paraSO.dc),...
        '_tau',num2str(paraSO.tau0*(paraSO.dtau^(paraSO.taut-1))),'_ISSymme',num2str(paraSO.IsSymme)];
elseif(Ni==1)
    paraSO=varargin{1};
end

%% Initalize state and define Hamiltonian
[Hx,Hy,Hz]=HamiltonianKitaHesen(paraSO.Jh,paraSO.Jk,paraSO.hx,paraSO.hz);
A{1}=InitialPEPS(paraSO.dc,d); A{2}=A{1};
LM{1,1}=ones(paraSO.dc,1)/sqrt(paraSO.dc); LM{2,1}=LM{1}; LM{3,1}=LM{1};
I4=ones(d*d,1);
paraSO.tau=paraSO.tau0*(paraSO.dtau.^(0:paraSO.taut-1));
Ntau=paraSO.taut;
Eb=zeros(paraSO.time,3); Etot=zeros(paraSO.time,1);
OBsimp.Ez=zeros(Ntau,1); OBsimp.Ex=OBsimp.Ez; OBsimp.Ey=OBsimp.Ez;
OBsimp.Ef=zeros(Ntau,1);

%%Iteration over tau
for nt=1:Ntau
    fprintf('Simulation with tau = %g  \n',paraSO.tau(nt));
    %% Prepare the evolution gates
    [ULx,URx]=ExpSvdU(Hx,paraSO.tau(nt));
    [ULy,URy]=ExpSvdU(Hy,paraSO.tau(nt));
    [ULz,URz]=ExpSvdU(Hz,paraSO.tau(nt));
    tob=0;
    for t=1:paraSO.time
        %% Evolution
        [A{1}]=EvolveTensor(A{1},ULx,2); [A{2}]=EvolveTensor(A{2},URx,2);
        LM{1}=kron(LM{1},I4);
        [A,LM]=SuperOrthogonalizationHoneycomPEPS(A,LM,paraSO.dc);
        
        [A{1}]=EvolveTensor(A{1},ULy,3); [A{2}]=EvolveTensor(A{2},URy,3);
        LM{2}=kron(LM{2},I4);
        [A,LM]=SuperOrthogonalizationHoneycomPEPS(A,LM,paraSO.dc);
        
        [A{1}]=EvolveTensor(A{1},ULz,4); [A{2}]=EvolveTensor(A{2},URz,4);
        LM{3}=kron(LM{3},I4);
        [A,LM]=SuperOrthogonalizationHoneycomPEPS(A,LM,paraSO.dc);
        
        %% Observation
        IsOb=(rem(t,paraSO.dtime)==0);
        if(IsOb)
            tob=tob+1;
            [Eb(tob,:)]=Observation(A,LM,{Hx,Hy,Hz});
            Etot(tob)=sum(Eb(tob,:))/2;
        end
        
        IsBreak=(IsOb && t>paraSO.time0 && tob>1.1 && abs(Etot(tob)-Etot(tob-1))<paraSO.Eps);
        if(IsBreak)
            fprintf('Converged at the %g-th iteraction, Etot = %g \n',t,Etot(tob))
            break;
        elseif(rem(t,paraSO.dt_print)==0)
            fprintf('At the %g-th iteration, Etot = %g \n',t,Etot(tob))
        end
    end
    OBsimp.Ex(nt)=Eb(tob,1);
    OBsimp.Ey(nt)=Eb(tob,2);
    OBsimp.Ez(nt)=Eb(tob,3);
    OBsimp.Ef(nt)=Etot(tob);
end

save([paraSO.PATHtree,paraSO.EXP,'.mat'],'paraSO','A','LM','ULx','URx','ULy','URy','ULz','URz','Hx','Hy','Hz','OBsimp');
end

function A=InitialPEPS(chi,d)
A=randn(d,chi,chi,chi);
A=A+permute(A,[1,3,4,2])+permute(A,[1,4,2,3]);
end

function [UL,UR]=ExpSvdU(H,tau)
d=2; D=d*d;
% U=reshape( permute( reshape(expm(-tau*H),[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );
U=reshape( permute( reshape(eye(size(H))-tau*H,[d,d,d,d]),[1,3,2,4] ),[d*d,d*d] );

[UL,LM,UR]=svd(U);
UL=permute( reshape( UL*sqrt(LM),[d,d,D] ),[1,3,2] ); 
UR=permute( reshape( conj(UR)*sqrt(LM),[d,d,D] ),[1,3,2] );
end

function A=SymmetrizeA(A) %#ok<DEFNU>
A=(A+permute(A,[1,2,5,4,3,6]))/2;
A=(A+permute(A,[1,2,3,6,5,4]))/2;
end

function rho2=ObtainRho2(A1,A2,invLMx,x)
S=size(A1);
A1=reshape( permute(A1,[2:x-1,x+1:4,1,x]),[prod(S([2:x-1,x+1:4])),S(1)*S(x)] );
A1=reshape( permute( reshape( A1'*A1,[S(1),S(x),S(1),S(x)] ),[1,3,2,4] ),[S(1)*S(1),S(x)*S(x)] );
A2=reshape( permute(A2,[2:x-1,x+1:4,1,x]),[prod(S([2:x-1,x+1:4])),S(1)*S(x)] );
A2=reshape( permute( reshape( A2'*A2,[S(1),S(x),S(1),S(x)] ),[2,4,1,3] ),[S(x)*S(x),S(1)*S(1)] );
rho2=reshape( permute( reshape( A1*diag(kron(invLMx,invLMx))*A2,ones(1,4)*S(1) ),[1,3,2,4] ),[S(1)*S(1),S(1)*S(1)] );
rho2=rho2/trace(rho2);
end

function [Eb]=Observation(A,LM,H)
Eb=zeros(1,3);
[A1]=Obsorb_lambda_PEPS(A{1},LM);
[A2]=Obsorb_lambda_PEPS(A{2},LM);
for x=2:1:4
    rho2=ObtainRho2(A1,A2,LM{x-1}.^(-1),x);
    Eb(x-1)=trace(H{x-1}*rho2);
end
end

function [T]=EvolveTensor(T,U,x)
% Evolve on the x-th bond (not x-th virtual bond)
S=size(T); [d,D,~]=size(U);

T=reshape( reshape(U,[d*D,d]) * reshape(T,[S(1),prod(S(2:4))]),[d,D,S(2:4)] );
T=permute(T,[1,3:x,2,x+1,x+2:5]);
S(x)=S(x)*D;
T=reshape(T,S);
end
